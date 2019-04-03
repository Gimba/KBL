#! /usr/bin/env python

# a nuisance is that crystal structures of the same protein use a different numbering for residues (e.g. C2-domain from 1d7p and C2-domain taken from the complex 1iqd). This script adjusts the numbering in this case.

import os, re

for file_name in os.listdir("."):
    if not 'histo' in file_name:
        number = re.findall(r'\d+', file_name)
        if number:
            new_number = int(number[-1]) - 3
            if new_number < 1:
                os.remove(file_name)
                continue
            new_file_name = file_name.replace(number[-1], str(new_number))
            os.rename(file_name, new_file_name)
