"""
TODO: Please Give me a name
======================
And a brief description.

Input files:

 - DIFF: file with a set of instruction which area applied to the dna sequences in
         the REFERENCE file.
 - REFERENCE:  File with reference DNA sequences.

Output files:

 - DIFF_OUT: DIFF file with reference dna after the instructions are applied.
 - LOG: Log file.  This includes all instances where there was an instruction 
        that could not be handled properly.
 
Handling instructions.  Reverse sort instructions by index and apply the
instruction backwards.  Assume that two instructions do not affect the same range
of characters (there is no reason for this I can think of).

Assume:

 - No replacing multiple bases.
 - 'd' means delete
 - Convert the bases to upper case for consistency

Output file
-------------
The fixed DNA string is saved as a new field at the end of the ``DIFF.tsv`` file.

So the fixed file looks like this::

  LOC;          SVTYPE;    DIFF;                    fixed
  Chr1:73857	ALU	       t126a,c128a,d134-282;     GGCCGGGCGCGGTGGC...

Things to check
----------------
Take a close look at insert behavior.  Are we off by 1.

License
----------
This code is published under the MIT license.

Copyright (c) 2017 Chip Webber
Copyright (c) 2017 Daniel Webber

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import os
import os.path as osp
import sys
from string import upper
import string
import re
from copy import copy
import pdb

DIFF = 'DIFF.tsv'
DIFF_OUT = 'DIFF.fixed.tsv'
REFERENCE = 'Reference.tsv'
LOG = 'log.txt'

# Insert pattern
P_INSERT = re.compile('i(\d+)([actgn]+)')
P_REPLACE = re.compile('([actg])(\d+)([actgn])')
P_N_ONE = re.compile('[n](\d+)')
FL = open(LOG, 'w')

with open(REFERENCE, 'r') as f:
    lines = f.readlines()

# len(ALU) = 312
# len(LINE1) = 5403
ALUJ = lines[0].split('\t')[1].strip()
ALUS = lines[1].split('\t')[1].strip()
ALUY = lines[2].split('\t')[1].strip()
LINE1 = lines[3].split('\t')[1].strip()

def check_base(base):
    base = list(base)
    for b in base:
        if string.upper(b) not in ['A', 'C', 'T', 'G', 'N']:
            raise ValueError("Invalid Base")


def parse_instr_set(instr_set, loc):
    instrs = [parse_instr(instr, loc) for instr in instr_set.split(',')]
    # Drop bad instructions.  Bad instructions are None, good instructions
    # are dicts
    valid_instrs = []
    for instr in instrs:
        if instr:
            valid_instrs.append(instr)

    # Handle the trivial case of no changes, or only one change
    if len(valid_instrs) == 1:
        return valid_instrs
    # Reverse sort by idx1
    sorted_instr = sorted(valid_instrs, key=lambda k: k['idx1'], reverse=True)
    return sorted_instr

def parse_instr(instr, loc):
    """
    Parse one instruction.

    Params:
      instr (string): a single instruction such as ``c174g`` or ``d134-282``
      loc (string): LOC field from DIFF.tsv, included for debug traces.
    """
    # Default values in dict
    idx1 = None
    idx2 = None
    base_old = None
    base_new = None

    instr = instr.strip()
    if 'NoDifference' in instr:
        itype = 'skip'
    elif len(instr) == 0:
        FL.write("Empty (Invalid) instruction at: %s\n" % loc)
        return None
    elif instr[0] == 'd' and '-' not in instr:
        itype = 'delete_one'
    elif instr[0] == 'd' and '-' in instr:
        itype = 'delete_multiple'
    elif instr[0] == 'n' and '-' not in instr:
        itype = 'n_one'
    elif instr[0] == 'n' and '-' in instr:
        itype = 'n_multiple'
    elif instr[0] == 'i':
        itype = 'insert'
    elif instr[0] in ['a', 'c', 'g', 't']:
        itype = 'replace'
    else:
        msg = "loc: %s  Invalid instruction: %s\n" % (loc, instr)
        FL.write(msg)
        #raise ValueError(msg)
        return None

    # Now find the index (for sorting) of the instruction
    if itype == 'skip':
        pass
    elif itype == 'delete_one':
        idx1 = int(instr[1:]) - 1
    elif itype == 'delete_multiple':
        idx1, idx2 = instr[1:].split('-')
        idx1 = int(idx1) - 1
        idx2 = int(idx2) - 1
    elif itype == 'insert':
        # example: i252gcagtcc
        m = P_INSERT.match(instr)
        if m is None:
            pdb.set_trace()
            FL.write("loc: %s Invalid insert instruction: %s\n" % (loc, instr))
            return None
            #raise ValueError("Invalid Insert statement: %s" % instr)
        idx1 = int(m.groups()[0]) - 1
        base_new = upper(m.groups()[1])
        check_base(base_new)
    elif itype == 'replace':
        m = P_REPLACE.match(instr)
        if m is None:
            FL.write("loc: %s Invalid replace instruction: %s\n" % (loc, instr))
            return None
        base_old = upper(m.groups()[0])
        idx1 = int(m.groups()[1]) - 1
        base_new = upper(m.groups()[2])
        for base in [base_old, base_new]:
            check_base(base)
    elif itype == 'n_one':
        #idx1 = int(instr[1:]) - 1
        m = P_N_ONE.match(instr)
        if m is None:
            FL.write("loc: %s Invalid n_one instruction: %s" % (loc, instr))
            return None
        idx1 = int(m.groups()[0]) - 1
    elif itype == 'n_multiple':
        idx1, idx2 = instr[1:].split('-')
        idx1 = int(idx1) - 1
        idx2 = int(idx2) - 1
    else:
        msg = "Invalid instruction type: %s" % instr
        raise ValueError(msg)

    dd = {'itype': itype,
          'idx1': idx1,
          'idx2': idx2,
          'base_old': base_old,
          'base_new': base_new,
          'raw': instr}
    return dd.copy()

def apply_instr(instr, dna, loc):
    """
    Apply the instruction to the dna sequence.

    Don't work with strings in python.  Work with lists then convert them
    to strings when you need them.  Strings are immutable and don't have all
    the cool methods lists have.

    Parameters:
      instr (dictionary): instruction dictionary
      dna (str): DNA string cleaned up from Reference.tsv
      loc (str): LOC field from DIFF.tsv for debug

    Returns:
      str: DNA string with instruction applied
    """
    dd = instr
    dna = list(dna)
    orig_dna = copy(dna)
    idx1 = dd['idx1']
    idx2 = dd['idx2']
    if dd['base_new']:
        base_new = dd['base_new']
        base_new = string.upper(base_new)
        base_new = list(base_new)
    if dd['itype'] == 'skip':
        return dna
    elif dd['itype'] == 'replace':
        if len(dna) < idx1:
            msg = "loc: %s  idx1 out of range in: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        if dna[idx1] != dd['base_old']:
            #pdb.set_trace()
            #raise ValueError("Invalid Replace Instruction")
            # TODO: figure out logging
            msg = "loc: %s  Invalid replace instruction: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        dna[idx1:idx1] = base_new
    elif dd['itype'] == 'delete_one':
        if len(dna) < idx1:
            msg = "loc: %s  idx1 out of range in: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        dna.pop(idx1)
    elif dd['itype'] == 'delete_multiple':
        if len(dna) < idx1:
            msg = "loc: %s  idx1 out of range in: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        if len(dna) < idx2:
            msg = "loc: %s  idx2 out of range in: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        dna[idx1:idx2] = []
    elif dd['itype'] == 'insert':
        # list.insert() inserts at a position.
        # we want to insert after a position
        if len(dna) < idx1:
            msg = "loc: %s  idx1 out of range in: %s\n" % (loc, instr['raw'])
            FL.write(msg)
            return dna
        dna[idx1+1 : idx1+1] = base_new
    elif dd['itype'] == 'n_one':
        dna[idx1] == 'n'
    elif dd['itype'] == 'n_multiple':
        orig_len = len(dna)
        n_len = idx2+1 - idx1
        assert n_len > 0
        n_list = ['n'] * n_len
        dna[idx1:idx2+1] = ['n'] * n_len
        new_len = len(dna)
        assert orig_len == new_len

    dna = "".join(dna)
    return dna

def fix_line(line):
    """
    Apply all the fixes to one line from the DIFF.tsv file
    """
    line = line.strip()
    if line.startswith('#LOC'):
        return "#LOC\tSVTYPE\tDIFF\tFIXED"
    loc, svtype, instr_orig = line.split('\t')
    instructions = parse_instr_set(instr_orig, loc)

    if upper(svtype) == 'ALUJ':
        dna = ALUJ
    elif upper(svtype) == 'ALUS':
        dna = ALUS
    elif upper(svtype) == 'ALUY':
        dna = ALUY
    elif upper(svtype) == 'LINE1':
        dna = LINE1
    else:
        raise ValueError("Invalid svtyp:%s" % svtype)

    for instr in instructions:
        dna = apply_instr(instr, dna, "%s %s" % (loc, svtype))

    line = "%s\t%s" % (line, dna)
    return line

with open(DIFF, 'r') as f_in:
    with open(DIFF_OUT, 'w') as f_out:
        for line in f_in:
            new_line = fix_line(line)
            f_out.write("%s\n" % new_line)
FL.close()
