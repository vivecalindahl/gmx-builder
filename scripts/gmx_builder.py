#! /usr/bin/env python2.7

# ======================================================
# Imports
# ======================================================


# Execute global startup script
import os

startup = os.environ.get("PYTHONSTARTUP")
if startup and os.path.isfile(startup):
    execfile(startup)

import argparse
import subprocess
from random import randint

# ======================================================
# Executing in the shell
# ======================================================

def run_in_shell(command):
    if command.count('|') >= 1:
        return run_in_shell_piped(command)

    proc = subprocess.Popen((command).split(), preexec_fn=os.setsid, close_fds = True,
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout, stderr  = stdout.strip(), stderr.strip()

    error = proc.returncode != 0
    if error:
        print stdout, stderr
        sys.exit("Failed to execute command: " + command + 
                 " (in: " + os.getcwd() + ")")

    return stdout

def run_in_shell_piped(command):
    commands = command.split('|') # this could fail

    # Chain piped commands together
    prev_proc = None;
    proc_piped=[]

    for c in commands:
        prev_stdout = None
        if prev_proc != None:
            prev_stdout = prev_proc.stdout

        proc = subprocess.Popen(c.split(), preexec_fn=os.setsid, close_fds = True,
                                           stdout = subprocess.PIPE, stderr= subprocess.PIPE,
                                           stdin = prev_stdout)
        proc_piped.append(proc)
        prev_proc = proc
    stdout, stderr = proc_piped[-1].communicate()
    stdout, stderr = stdout.strip(), stderr.strip()

    error = proc.returncode != 0
    if error:
        print stdout, stderr
        sys.exit("Failed to execute command: " + command)

    return stdout

def run_in_shell_allow_error(command):
    proc = subprocess.Popen((command).split(), preexec_fn=os.setsid, close_fds = True,
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.strip()

    error = proc.returncode != 0

    return stdout, stderr, error

# ======================================================
# Files and paths
# ======================================================

def absolute_path(path):
    stdout = run_in_shell('readlink -f ' + path)
    return stdout

def clone_directory(srcdir, outdir):
    contents=' '.join(['/'.join([srcdir,item]) for item in os.listdir(srcdir)])
    run_in_shell('mkdir -p ' + outdir)
    run_in_shell(' '.join(['cp -r', contents, outdir]))

def remove_temporary_files():
    tmp = ' '.join([ f for f in os.listdir(os.getcwd())  if f.startswith('#') or f.endswith('~') or f.startswith('_tmp') ])
    run_in_shell('rm -rf ' + tmp)

def concatenate_files(filenames, outfilename):
    with open(outfilename, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def read_lines(inputfile):
    with open(inputfile) as f:
        lines=f.read().splitlines()
    return lines

def write_lines(lines, outputfile):
    # will overwrite existing file
    with open(outputfile, 'w') as f:
        for line in lines:
            f.write(line + '\n')

def check_path_is_clean(path, forceful=False):
    if os.path.exists(path) and not forceful:
        sys.exit(path + ' already exists. Use -f to force overwrite.')
    else:
        run_in_shell('rm -rf ' + path)        

# ======================================================
# gmx related
# ======================================================

# gro-file format info
gro_header_size = 2
gro_footer_size = 1

# The keys and their order
gro_keys = ['resid', 'residue_name', 'atom_name', 'atomid', 
            'x', 'y', 'z',
            'vx', 'vy', 'vz']

# The number of positions for the value of each key 
gro_npositions = {'resid':5, 'residue_name':5, 'atom_name':5, 'atomid':5,
                  'x':8, 'y':8, 'z':8,
                  'vx':8, 'vy':8, 'vz':8}
gro_table_columns = {k:v  for k, v in  zip(gro_keys, range(len(gro_keys)))}

def gro_table_column(key):

    if key not in set(gro_keys):
        raise ValueError("results: status must be one of %r." % gro_keys)

    return gro_table_columns[key]

# Return string table with the gro file data, headers excluded
def read_gro_table(gro):
    # Read gro-file as a list of strings
    lines = read_lines(gro);

    # Remove data headers
    lines = lines[gro_header_size:-gro_footer_size]

    table=[]
    for line in lines:
        row = []
        position = 0
        for k in gro_keys:
            npos = gro_npositions[k]
            value = line[position:position+npos].strip()
            position += npos
            row.append(value)
            
        table.append(row)

    return table

def read_gro_coordinates(gro):
    table = read_gro_table(gro)

    # Extract coordinates as numpy array
    cols = [gro_table_columns[coord] for coord  in ['x', 'y', 'z']]
    coords=[]
    for row in table:
        xyz = []
        for col in cols:
            xyz.append(float(row[col]))
        coords.append(xyz)

    coords = np.array(coords)

    return coords

def merge_dicts(dlist):
    merged={}
    for d in dlist:
        # An dictionary further back in the list has precedence over earlier ones
        for k, v in d.iteritems():
            merged[k] = v

    return merged

def merge_mdps(mdplist):
    return merge_dicts(mdplist)

def write_mdp(mdp, outputfile):
    lines =  [ ' = '.join([k,v]) for k, v in  mdp.iteritems() ]

    # Here just sort alphabetically, something else is surely optimal...
    # Could also tweak the output format.
    lines.sort()

    write_lines(lines, outputfile)

def write_selections(sel, outputfile):
    # selections are a list of strings. Statements should be separated by ';'.
    lines = [s + ';' for s in sel if not s.endswith(';') ]
    write_lines(lines, outputfile)

def make_tpr(templatedir, nomdp=False, indexfile=False, out='topol.tpr'):
    # Work in a temporary directory.
    startdir=absolute_path(os.getcwd())
    templatedir=absolute_path(templatedir)
    workdir='./_tmpdir' + str(randint(1,1e4))
    clone_directory(templatedir, workdir)
    os.chdir(workdir)

    if nomdp:
        run_in_shell('touch grompp.mdp')

    extra_args=''   
    if indexfile:
        extra_args += ' -n'

    # Allow warnings but print if there were some.
    grompp = ' '.join(['gmx grompp', extra_args])
    stdout, stderr, error = run_in_shell_allow_error(grompp)
    if error:
        extra_args += ' -maxwarn 10'
        grompp = ' '.join(['gmx grompp', extra_args])
        stdout = run_in_shell(grompp)
        print("'gmx grompp' generated warnings for " + templatedir)

    # Only keep the tpr-file
    run_in_shell('mv topol.tpr ' + startdir + '/' + out)
    os.chdir(startdir)
    run_in_shell('rm -r ' + workdir)    

def make_index_file_from_selections(sel_file, tpr='topol.tpr'):
    # If grompp gets a custom index file it will be unaware about the 
    # otherwise generated default groups, e.g. "system", which are
    # generally needed e.g. to specify the  temperature coupling groups.
    # So, generate a default index file first, then add the actual selections.

    defaults='defaults.ndx'

    run_in_shell_piped('echo q | gmx make_ndx -o ' + defaults + ' -f ' +  tpr)
    # Note: the shell commands can be a bit sensitive. If q --> "q" above, it fails.
    # It corresponds to writing 'echo \"q\" | [...]' directly in the terminal.

    selected='selected.ndx'
    run_in_shell(' '.join(['gmx select','-s', tpr, '-sf', sel_file, '-on', selected]))

    concatenate_files([defaults, selected], 'index.ndx')

    run_in_shell(' '.join(['rm', defaults, selected]))

def make_run_template(setup, mdp_settings, outdir,  selections=[]):
    clone_directory(setup, outdir)

    mdp_out = outdir + '/grompp.mdp'
    write_mdp(mdp_settings, mdp_out)

    startdir = os.getcwd()
    os.chdir(outdir)

    # Generate an index file from the selections, if given. 
    have_selections = len(selections) > 0
    if have_selections:

        sel_out = '_tmp-sel-' + str(randint(1,1e4)) + '.txt'
        write_selections(selections, sel_out)

        # a tpr-file is assumed to be already available in the setup directory
        make_index_file_from_selections(sel_out)
        remove_temporary_files()

    # The ultimate test: is it possible to make a tpr from this?
    make_tpr(outdir, indexfile = have_selections)
    os.chdir(startdir)
