#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 13:44:35 2020

@author: ceballos
"""

from subprocess import check_call, STDOUT
import shlex

def run_comm(comm, msg=""):
    print(msg)
    print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except Exception as mess:
        print(mess)
        raise