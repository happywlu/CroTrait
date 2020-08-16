# -*- coding: utf-8 -*-
# @Author: Lenovo1
# @Date:   2020-08-14 18:07:36
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-14 23:46:35

import os
import sys
import subprocess

star_dir = "example_sequence"

ap = os.path.abspath("%s" % star_dir)

ab = subprocess.call(["where", "blastn"])

print(ab)

ap1 = os.path.abspath("%s" % "Cronobacter_OACs.fasta")
m = os.path.split(ap)

print(ap1)

print(os.path.split("Cronobacter_OACs.fasta"))