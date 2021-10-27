#!/usr/bin/python
#
# python send_email [RECEIVER_EMAIL]
#
# Send email after all study checks complete...
#

# import necessary libraries
import json
import glob
import re
import os
from os import path
import subprocess
from subprocess import call
from pathlib import Path
import warnings


#get input arguments
import sys, getopt

# set up email notification
import smtplib, ssl

port = 465  # For SSL
smtp_server = "smtp.gmail.com"
sender_email = "noreply.incdata@gmail.com"  # Enter your address
password = "XXXXXX".  # password hashed for public view

receiver_email=sys.argv[1]

f=open("email.txt","r")
if f.mode=="r":
	message = f.read()
f.close()


context = ssl.create_default_context()
with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
    server.login(sender_email, password)
    server.sendmail(sender_email, receiver_email, message)
