#! usr/bin/env python

# ## PIPELINE: inc_dcm2bids.py
# ## USAGE: python3 inc_dcm2bids --template=<templatefile> --trange=<days> [OPTIONS]
#    * requires python3, freesurfer, FSL (calls FSL via python subprocess)
#
# ## Author(s)
#
# * Amy K. Hegarty, Intermountain Neuroimaging Consortium, University of Colorado Boulder
# * University of Colorado Boulder
#
# ## Product
#
# FSL Pipelines
#
# ## License
#
# <!-- References -->
# [FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki
# [pybids]: Yarkoni et al., (2019). PyBIDS: Python tools for BIDS datasets. Journal of Open Source Software, 4(40), 1294, https://doi.org/10.21105/joss.01294
#           Yarkoni, Tal, Markiewicz, Christopher J., de la Vega, Alejandro, Gorgolewski, Krzysztof J., Halchenko, Yaroslav O., Salo, Taylor, ? Blair, Ross. (2019, August 8). bids-standard/pybids: 0.9.3 (Version 0.9.3). Zenodo. http://doi.org/10.5281/zenodo.3363985
#
import os, sys, getopt, glob, re, bids, json, warnings, subprocess, multiprocessing
from subprocess import PIPE
from os import path
import pandas as pd
import os.path, time
import datetime as dt
from pathlib import Path
import collections
import smtplib, ssl

# ------------------------------------------------------------------------------
#  Show usage information for this script
# ------------------------------------------------------------------------------

def print_help():
  print("""
    Intermountain Neuroimaging Consortium Dicom 2 BIDS conversion
        Usage: """ + """ --template=<bids-template> --trange=<days> [OPTIONS]
        OPTIONS
          -h --help                   show this usage information and exit
          -i --template               [JSON] study template describing naming convention
          -t --trange                 time range (days) for dicom retreival (range starts today)
          -c --nifti-convert          (default false) run conversion to nifti for
                                      scanner dicoms. Output in bids format
          -l --subject-id-key         provide a key for M803 number and subject ids (overrides accession #)
                                      format: (csv) column 1: run(1)/skip(0)
                                                    column 2: dicom path (USRI/Study date and time
                                                    column 3: subject id
                                                    column 4: (optional) session id
          -w --suppress-warnings      (default false) suppress warings on terminal
    ** OpenMP used for parellelized execution of XXX. Multiple cores (CPUs) 
       are recommended (XX cpus for each fmri scan).
       
    ** see github repository for more information and to report issues: 
       https://github.com/intermountainneuroimaging/dcm2bids.git
        
        """)


# ------------------------------------------------------------------------------
#  Parse arguements for this script
# ------------------------------------------------------------------------------

def parse_arguments(argv):

    #intialize arguements
    print("\nParsing User Inputs...")
    runconvert = False
    warns = False
    pidkey = ""
    wd=os.popen('echo $HOME/scratch').read().rstrip()

    try:
      opts, args = getopt.getopt(argv,"hi:t:l:cw",["template=","trange=","help","nifti-convert","subject-id-key=","suppress-warnings"])
    except getopt.GetoptError:
      print_help()
      sys.exit(2)
    for opt, arg in opts:
      if opt in ("-h", "--help"):
         print_help()
         sys.exit()
      elif opt in ("-i", "--template"):
         bidstemplate = arg
         if not os.path.exists(bidstemplate):
           raise Exception("BIDS study template does not exist")
      elif opt in ("-t", "--trange"):
         trange = arg
      elif opt in ("-c","--nifti-convert"):
         runconvert = True 
      elif opt in ("-l","--subject-id-key"):  
         pidkey = arg 
      elif opt in ("-w","--suppress-warnings"):
         warns = True                           
    if 'bidstemplate' not in locals():
      print_help()
      raise Exception("Missing required argument -i [--template]")
      sys.exit()
    if 'trange' not in locals():
      print_help()
      raise Exception("Missing required argument -t [--trange]")
      sys.exit()

    class args:
      def __init__(self, wd, bidstemplate, trange, runconvert,warns,pidkey):
        self.wd = wd
        self.bidstemplate = bidstemplate
        self.trange = trange
        self.runconvert = runconvert
        self.suppresswarnings = warns
        self.pidkey = pidkey
        self.templates = "/projects/amhe4269/banichlab_ldrc_preproc/inc_resources/scanner_check/v2.0/"

    entry = args(wd, bidstemplate, trange, runconvert,warns,pidkey)

    return entry

# # ------------------------------------------------------------------------------
# #  Main Pipeline Starts Here...
# # ------------------------------------------------------------------------------

def worker(name,cmdfile):
    """Executes the bash script"""
    process = subprocess.Popen(cmdfile.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output, error = process.communicate()
    print(output)
    print(error)
    print("worker: " + name + " finished")
    return


# define functions
def list_diff(list1, list2): 
  return (list(set(list1) - set(list2))) 

def last_4chars(x):
  return(x[-4:])

# load study heuristic file...
def heuristic(entry):

  # load study template
  with open(entry.bidstemplate) as f:
    data = json.load(f)

  # get study specifics from heuristic file
  heuristic.studyname=data['Acquisition'][0]['Study'][0]['name']
  heuristic.scannerpath=data['Acquisition'][0]['Study'][0]['scanner_regexp']

  # check bids path if defined
  if 'bids' in data['Acquisition'][0]['Study'][0]:
    heuristic.bidspath=data['Acquisition'][0]['Study'][0]['bids']
    if heuristic.bidspath != "":
      #check bids directory exists...
      if not path.exists(heuristic.bidspath):
        raise Exception([heuristic.bidspath + " does not exist!"] )
  else:
    bidspath=[]

  # get regular expressions for T1w for ascension command
  heuristic.subregexp=data['Acquisition'][0]['Subject'][0]['scanner_regexp']
  heuristic.sesregexp=data['Acquisition'][0]['Session'][0]['scanner_regexp']
  heuristic.t1wregexp=data['Acquisition'][0]['anat'][0]['scanner_regexp']

  # check scanner path exists...
  if not path.exists(heuristic.scannerpath):
    raise Exception([heuristic.scannerpath + " does not exist!"] )

  heuristic.strexp=heuristic.scannerpath + '/' + heuristic.subregexp + '/' + heuristic.sesregexp + '/' + heuristic.t1wregexp + '_????/*0001-1.dcm'
  heuristic.data=data

def get_pidkey(entry):

  if entry.pidkey:
    key=pd.read_csv(entry.pidkey, header=None)
    print("Using Key:")
    print(key)
    key.drop(key[key[0] == 0].index, inplace = True)
  else:
    key=pd.DataFrame()  

  return key


# Get Subject / Session Set
def new_scans(entry):

  # check time ids for all recent file transfers - if within time domain report...
  today = dt.datetime.now().date()
  start_date = today - dt.timedelta(days=int(entry.trange))
  print("Searching for DICOMS: " + start_date.strftime("%Y-%m-%d")  + " to " + today.strftime("%Y-%m-%d")  + "... ")

  # pull all recent files...
  flag_newscan = False
  alltxt=""
  key = get_pidkey(entry)

  #add entry to log
  bidspath = heuristic.data['Acquisition'][0]['Study'][0]['bids']
  f = open(bidspath+"/tmp/dcm2bids.key", "a")
  f.write("Entry: "+today.strftime("%m/%d/%Y, %H:%M:%S")+"\n")
  f.close()

  files=subprocess.check_output("cd "+heuristic.scannerpath+"; find M803*/Study* -maxdepth 0 -mtime -"+entry.trange,shell=True)
  files = files.decode().strip().split('\n')

  if files[0]:
    for f in files:
      
      if 'layout' not in locals():
        layout = make_bidsobject(bidspath)

      # ppath = s.join(fpath)
      ppath=heuristic.scannerpath+"/"+f
      print("\nAcquisition...")
      print(ppath)
      print(" ")

      study = heuristic.data['Acquisition'][0]['Study'][0]['name']
      subdigit = heuristic.data['Acquisition'][0]['Subject'][0]['digits']
      sesdigit = heuristic.data['Acquisition'][0]['Session'][0]['digits']
      scannerID=ppath.split('/')[-2]
      s="/"
      scannertimedate=s.join(ppath.split('/')[-2:])

      filetime = dt.datetime.fromtimestamp(os.path.getmtime(ppath))
      scandate=filetime.date().strftime("%Y-%m-%d")

      if key.empty:
        filepath=ppath + '/' + heuristic.t1wregexp + '_????/*0001-1.dcm'
        
        dicomfile=subprocess.check_output("ls "+filepath,shell=True)
        dicomfile=dicomfile.decode().strip().split('\n')

        print("Pulling Subject ID from DICOM accession...")
        process = subprocess.Popen(['./accession.sh',dicomfile[0]], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out, err = process.communicate()
        out = out.strip('\n')
        print(err)

        out = out.replace('rray','')  # special rule for rray study...
        out = out.replace('cbd','')  # special rule for cbdx study...
        out = out.replace('cwb','')  # special rule for cwb study...

        if '/' in out:
          pid = out.split('/')[1].zfill(subdigit)      # pull subject ID from accession number
          if len(out.split('/')) > 2:
            ses = out.split('/')[2]                      # pull session ID from accession number
            ses = ses.strip('s').strip('S').zfill(sesdigit)
          else:
            ses = "none"
        else:
          pid = out
          ses = "none"

        if not pid:   # if pid is empty (issue with accession)
          pid="unknown"

      else:
        index = key.index[key[1] == scannertimedate].tolist()
        print(index)
        if not index:
          continue
        pid = str(key.loc[index[0],2])

        if len(key.columns) > 3:
          ses = key[3][index[0]]
        else:
          ses = "none"

      print("Running for sub-" + str(pid) + " ses-" + str(ses))
      ff = open(bidspath+"/tmp/dcm2bids.key", "a")
      ff.write("Convert sub-" + str(pid) + " ses-" + str(ses)+"\n")
      ff.close()

      # match image with template info
      new_aquisitions = make_bidsname(entry,ppath,pid,ses,heuristic.data,layout)

      # convert new aquisitions to nifit
      
      nifti_convert(entry,new_aquisitions,heuristic.data)

      # generate report text
      rptext = make_textreport(study,scannertimedate,pid,ses,scandate)
      print(rptext)

      # add to summary scanner email (for all studies)
      if not os.path.exists("email.txt"):
        fo=open("email.txt","w+")
        fo.write("Subject: " + "[inc_scanner_report] " + "Scanner Check : date range " + str(start_date) + " to " + str(today) + "\n\n" )
        fo.close() 

      fo = open('email.txt', 'a')
      fo.write(rptext + "\n\n")
      # Close the file
      fo.close()     

      alltxt=alltxt+rptext + "\n\n"

      flag_newscan = True

  if not flag_newscan:
    print("No aquisitions in selected time range")

  if flag_newscan:
    sendemail(alltxt,heuristic.data,entry)

  # END new_scans


def make_bidsobject(bidspath):

  if not os.path.exists(bidspath) or not os.path.exists(bidspath +'/' + 'dataset_description.json'):
    os.makedirs(bidspath,exist_ok=True)
    
    # make dataset_description file...
    
    data = {
      'Name': 'Intermountain Neuroimaging Consortium Dataset',
      "BIDSVersion": "1.1.1",
      "CodeURL": "",
      "HowToAcknowledge": ""}

    with open(bidspath + '/' + 'dataset_description.json', 'w') as outfile:
        json.dump(data, outfile, indent=2)

  print("\n...Loading BIDS Object \n")

  bids.config.set_option('extension_initial_dot', True)

  layout = bids.BIDSLayout(bidspath, derivatives=False, absolute_paths=True)

  return layout

  # END make_bidsobject

def make_bidsname(entry,filepath,pid,ses,data,layout):

  # loop through all possible scanner images to convert (e.g. anat, func, dwi, ...)
  imgtype = ['anat', 'func', 'dwi', 'fmap']
  acq_list = [['dicomdir'],['bidsfile']]

  for t in imgtype:

    # all images of type t
    if t not in data['Acquisition'][0]:
      continue

    modalimgs=data['Acquisition'][0][t]

    # all images of imgtype "t" (e.g. loop through all func images)
    for img in modalimgs:

      name = img['name']
      scanner_regexp = img['scanner_regexp']
      input_regexp = filepath + '/' + scanner_regexp + '_????'

      directories = glob.glob(input_regexp)

      if 'nAcquisitions' in img:
        nruns = int(img['nAcquisitions'])
      else:
        nruns = 1
      
      # error checking here
      dcm_errorcheck(directories,img,t,nruns)   # check for errors in dicom set based on study template
                                              # return misstext, dupltext, incompltext for report...


      # Fix any issues...
      if len(directories) > nruns:
        directories.sort(key=os.path.getctime)
        directories = directories[-nruns:]       # if too many runs are logged, use the most recently collected set
        # raise Exception("Number of matching aquisitions exceeds expected scans")
        

      elif len(directories) < nruns:
        if not entry.suppresswarnings:
          print("**** \n\nNumber of matching aquisitions is less than expected scans:\n" + input_regexp + "\n Expected: "+ str(nruns) + " ... Found: " + str(len(directories)) + "\n\n***")
          warnings.warn("**** \n\nNumber of matching aquisitions is less than expected scans\n\n***")

      directories = sorted(directories, key = last_4chars) 

      # get match dicom name to bids name
      r=1
      for inputdir in directories:

        # Define the pattern to build out bids file structure
        pattern = data['Acquisition'][0]['Study'][0]['bids_pattern'] 
        # pattern = "sub-{subject}/[ses-{session}/][{type}/]sub-{subject}[_ses-{session}][_task-{task}][_acq-{acquisition}][_rec-{reconstruction}][_run-{run}][_echo-{echo}][_dir-{direction}][_space-{space}][_desc-{desc}]_{suffix}.nii.gz",

        # Populate values for pattern
        ent = img.copy()  # pull pattern values directly from json file... (e.g. task, aquisition, suffix)
        
        ent['subject'] = pid
        if ses != "none":
          ent['session'] = ses
        else:
          ent['session'] = []

        ent['type'] = t

        if "run" in img:
          if img["run"] == "n":
            ent['run'] = str(r).zfill(2)

        r = r + 1

        bidsfile = layout.build_path(ent, pattern, validate=False, absolute_paths=False)

        dicomdir = inputdir

        acq_list[0].append(dicomdir)  # all the dicom file paths
        acq_list[1].append(bidsfile)  # all the bids filenames

  return acq_list

  # END make_bidsname

def nifti_convert(entry,acq_list,data):
  
  warningstxt=""

  if not hasattr(nifti_convert,"warnings"):
    nifti_convert.warnings = "Warnings: "

  # check to make sure dicom list and bids list is the same length
  if len(acq_list[0]) != len(acq_list[1]):
    raise Exception("Issue with matching DICOM and BIDS format: Contact support team")
  
  studyname = data['Acquisition'][0]['Study'][0]['name']
  bidspath = data['Acquisition'][0]['Study'][0]['bids']
  r=1
  jobs=[];
  for i in range(1,len(acq_list[0])):
    dicomdir = acq_list[0][i]
    bidsfile = acq_list[1][i].replace(".nii.gz","")
    bidsname=bidsfile.split("/")[-1]
    # print(dicomdir + " to " + bidsfile)

    # check if output file already exists
    if os.path.exists(bidspath + '/' + bidsfile + ".nii.gz"):
      print("Warning: " + bidsfile + ".nii.gz" + " already exists")
      warningstxt=warningstxt + "\n" + bidsfile + ".nii.gz" + " already exists"
      continue

    # print list of aqusition + bids names
    s="/"
    dicomname=dicomdir.split("/")
    dicomname=s.join(dicomname[-3:])

    f = open(bidspath+"/tmp/dcm2bids.key", "a")
    f.write(dicomname + "\t" + bidsname +"\n")
    f.close()
    

    #check if nifti conversion should run
    if not entry.runconvert:
      continue

    # run dcm2niix on cluster...
    # call worker...

    # run dicom converter - pass output to bids directory
    subiden = [i for i in bidsname.split("_") if "sub" in i]
    sesiden = [i for i in bidsname.split("_") if "ses" in i]

    pid = subiden[0].split("-")[1]

    if sesiden:
      ses = sesiden[0].split("-")[1]
    else:
      ses=""

    imgdir=bidspath + "/tmp/data/bimages/" + pid + "/" + ses + "/Nifti"
    os.makedirs(imgdir, exist_ok=True)

    # prd_name=Path(dicomdir).stem[:-5]
    dicomname=dicomdir.split("/")
    prd_name=dicomname[-1]
    
    niftiimg = imgdir + '/' + prd_name
    
    print('Running:' + bidsfile + '.nii.gz')

    # add session paths
    cmdfile = imgdir + '/run_dcm2niix_' + str(r).zfill(2)

    cmd='cp ' + entry.templates + 'run_dcm2niix.sh ' + cmdfile
    subprocess.run(cmd.split())

    cmd='chmod u+x ' + cmdfile
    subprocess.run(cmd.split())

    sed1=('s^DICOM_PLACEHOLDER^' + dicomdir + '^')
    sed2=('s^BIDSFILE_PLACEHOLDER^' + bidspath + '/' + bidsfile + '^')
    sed3=('s^NIFTI_PLACEHOLDER^' + niftiimg + '^')
    subprocess.run(['sed','-i',sed1,cmdfile]) 
    subprocess.run(['sed','-i',sed2,cmdfile]) 
    subprocess.run(['sed','-i',sed3,cmdfile]) 


    # run script
    slurm_outfile=imgdir+"/dcm2niix.o%j"
    slurm_errfile=imgdir+"/dcm2niix.e%j"
    sbatchflags = "-q preemptable -p blanca-ccn -A blanca-ics-" + studyname + " -c 2 --job-name dcm2niix --wait --time=04:00:00 --mem=16G -o " + slurm_outfile + " -e " + slurm_errfile 
    cmd = 'sbatch ' + sbatchflags + ' ' + cmdfile

    name = "dcm2niix-" + str(r).zfill(2)
    p = multiprocessing.Process(target=worker, args=(name,cmd))
    jobs.append(p)
    p.start()

    r = r+1
    print(p)

    for job in jobs:
      job.join()  #wait for all distcorrepi commands to finish

  nifti_convert.warnings = nifti_convert.warnings + warningstxt

  # add intended for section in fieldmaps
  if entry.runconvert:
    intendedfor(acq_list,bidspath)

  # END nifti_convert

def intendedfor(acq_list,bidspath):
  
  # add intended for section in fmap json files...
  funcfiles = [i for i in acq_list[1] if "func/" in i]
  dwifiles =  [i for i in acq_list[1] if "dwi/" in i]
  fmapfiles = [i for i in acq_list[1] if "fmap/" in i]

  intendedfor_list = funcfiles + dwifiles
  sep='/'
  for i in range(0,len(intendedfor_list)):
  	file=intendedfor_list[i].split('/')
  	intendedfor_list[i]=sep.join(file[-2:])

  if fmapfiles:  # if convert list includes fmaps
    import time
    time.sleep(120) #sleep 30 seconds to make sure file shows up in new location

    for i in fmapfiles:

      i = i.replace(".nii.gz", ".json")
      with open(bidspath + '/' + i) as f:
          data = json.load(f)

      y={"IntendedFor":intendedfor_list}
      data.update(y)
      dat = collections.OrderedDict(data)
      dat.move_to_end('IntendedFor',last=True)

      data=dict(dat)

      with open(bidspath + '/' + i, 'w') as outfile:
          json.dump(data, outfile, indent=2)

  # END intendedfor


def dcm_errorcheck(directories,img,modality,nruns):

  # three primary error checks: missing, incomplete, duplicate runs

  duplimg = []
  missimg = []
  incompimg = []

  if not hasattr(dcm_errorcheck,"misstext"):
    dcm_errorcheck.misstext="Missing scans:"
  if not hasattr(dcm_errorcheck,"dupltext"):
    dcm_errorcheck.dupltext="Duplicate scans:"
  if not hasattr(dcm_errorcheck,"incomptext"):
    dcm_errorcheck.incomptext="Incomplete scans:"

  # check if directories matches number of runs
  if len(directories) > nruns:
    directories.sort(key=os.path.getctime)
    directories = directories[-nruns:]

    duplimg.append(modality + "-" + img['name'])

  elif len(directories) < nruns and len(directories) > 0:
    for rr in range(len(directories),nruns):
      missimg.append(modality + "-" + img['name'] + "_run-" + str(rr+1).zfill(2))
  elif len(directories) == 0:
    missimg.append(modality + "-" + img['name'])
  
  r=1;

  for inputdir in directories:

    # if length is defined...check this
    if 'length' in img:
      scan_length = img['length']

      # check scan length
      dcms = glob.glob(inputdir + "/*.dcm")
      if len(dcms) < scan_length:
        incompimg.append(modality + "-" + img['name'] + "_run-" + str(r).zfill(2))

    r=r+1

  nl="\n  "
  if missimg:
    dcm_errorcheck.misstext = dcm_errorcheck.misstext + "\n  " + nl.join(missimg)
  if duplimg:
    dcm_errorcheck.dupltext = dcm_errorcheck.dupltext + "\n  " + nl.join(duplimg)
  if incompimg:
    dcm_errorcheck.incomptext = dcm_errorcheck.incomptext + "\n  " + nl.join(incompimg)

  # END dcm_errorcheck


def make_textreport(study,scannerid,pid,ses,scandate):

  # store information for each subject ... run once per scan session

  # generate report text

  if study == "rray":
    if ses != "none":
  	  subjecttext="Subject: r" + str(pid) + "s" + str(ses)
    else:
      subjecttext="Subject: r" + str(pid)
  else:
  	subjecttext="Subject: " + str(pid) + " Session: " + str(ses)

  reporttxt = ("Study: " + study + "\n"
               "Scannerid: " + scannerid + "\n" + subjecttext + "\n"
               "Scan Date: " + scandate + "\n")
  if not hasattr(nifti_convert,"warnings"):
    nifti_convert.warnings = "Warnings:"
  if len(nifti_convert.warnings.split("\n")) < 2:
    nifti_convert.warnings = nifti_convert.warnings + " none"

  if len(dcm_errorcheck.misstext.split("\n")) < 2:
    dcm_errorcheck.misstext = dcm_errorcheck.misstext + " none"

  if len(dcm_errorcheck.incomptext.split("\n")) < 2:
    dcm_errorcheck.incomptext = dcm_errorcheck.incomptext + " none"

  if len(dcm_errorcheck.dupltext.split("\n")) < 2:
    dcm_errorcheck.dupltext = dcm_errorcheck.dupltext + " none"

  reporttxt=(reporttxt + nifti_convert.warnings + "\n" + dcm_errorcheck.misstext + "\n" + dcm_errorcheck.incomptext + "\n" + dcm_errorcheck.dupltext)

  delattr(nifti_convert,"warnings")
  delattr(dcm_errorcheck,"misstext")
  delattr(dcm_errorcheck,"incomptext")
  delattr(dcm_errorcheck,"dupltext")

  return reporttxt

# set up email notification
def sendemail(emailtxt,data,entry):
  

  port = 465  # For SSL
  smtp_server = "smtp.gmail.com"
  sender_email = "noreply.incdata@gmail.com"  # Enter your address
  password = "Buffs2021!"

  # get date info
  import os.path, time
  import datetime as dt
  today = dt.datetime.now().date()
  start_date = today - dt.timedelta(days=int(entry.trange))

  # # ... email text ....
  studyname = data['Acquisition'][0]['Study'][0]['name']
  message = "Subject: " + "[inc_scanner_report] " + studyname + ": date range " + str(start_date) + " to " + str(today) + "\n\n" + emailtxt

  for person in data['Study Contact']:
    receiver_email = person['email']

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, message)


def main(argv):
  
  # get user entry
  entry = parse_arguments(argv)

  os.makedirs(entry.wd, exist_ok=True)
  logdir = entry.wd + '/logs'

  # get participant bids path:
  heuristic(entry)

  # check for new scans
  new_scans(entry)


  # END main


if __name__ == "__main__":
    main(sys.argv[1:])
