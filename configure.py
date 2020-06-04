#! /usr/bin/env python
#-----------------------------------------------------------------------------------------
import argparse
import glob
import re

# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# Step 1. Prepare parser, add each of the arguments
parser = argparse.ArgumentParser()

user_directory = 'user/'
user_choices = glob.glob(user_directory + '*.F90')
user_choices = [choice[len(user_directory):-4] for choice in user_choices]

# system
parser.add_argument('-perseus',
                    action='store_true',
                    default=False,
                    help='Configure for `Perseus` cluster.')

parser.add_argument('-intel',
                    action='store_true',
                    default=False,
                    help='enable intel compiler')
parser.add_argument('-hdf5',
                    action='store_true',
                    default=False,
                    help='enable HDF5 & use h5pfc compiler')

parser.add_argument('-ifport',
                    action='store_true',
                    default=False,
                    help='enable IFPORT library (`mkdir` etc)')

mpi_group = parser.add_mutually_exclusive_group()
mpi_group.add_argument('-mpi',
                       action='store_true',
                       default=False,
                       help='enable mpi')
mpi_group.add_argument('-mpi08',
                       action='store_true',
                       default=False,
                       help='enable mpi_f08')

# user file
parser.add_argument('--user',
                    default=None,
                    choices=user_choices,
                    help='select user file')

# algorithms
parser.add_argument('--nghosts',
                    action='store',
                    default=5,
                    help='specify the # of ghost cells')

parser.add_argument('-absorb',
                    action='store_true',
                    default=False,
                    help='enable absorbing boundaries')

parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable DEBUG flag')

args = vars(parser.parse_args())

# Step 2. Set definitions and Makefile options based on above arguments

makefile_options = {}

makefile_options['USER_FILE'] = args['user']
makefile_options['USER_DIR'] = user_directory

makefile_options['COMPILER_COMMAND'] = ''
makefile_options['COMPILER_FLAGS'] = ''
makefile_options['PREPROCESSOR_FLAGS'] = ''

# specific cluster:
specific_cluster = False
if args['perseus']:
    specific_cluster = True
    args['intel'] = True
    args['mpi08'] = True
    args['mpi'] = False
    args['ifport'] = True
    makefile_options['COMPILER_FLAGS'] += '-xCORE-AVX2 '

# compilation command
if args['hdf5']:
    makefile_options['COMPILER_COMMAND'] += 'h5pfc '
    makefile_options['PREPROCESSOR_FLAGS'] += '-DHDF5 '
else:
    if ((not args['mpi']) and (not args['mpi08'])):
        makefile_options['COMPILER_COMMAND'] += 'gfortran '
    else:
        makefile_options['COMPILER_COMMAND'] += 'mpif90 '
if args['ifport']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DIFPORT '

# mpi version
if args['mpi']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI '
elif args['mpi08']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DMPI08 '

# debug
if args['debug'] and (not args['intel']):
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG -fcheck=all -fimplicit-none -fbacktrace '
if args['debug'] and args['intel']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DDEBUG '
    makefile_options['COMPILER_FLAGS'] += '-traceback '

# compilar (+ vectorization etc)
if args['intel']:
    makefile_options['MODULE'] = '-module '
    makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -xHost -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto '
else:
    makefile_options['MODULE'] = '-J '
    makefile_options['COMPILER_FLAGS'] += '-O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5 '

makefile_options['EXE_NAME'] = 'tristan-ff'

if args['absorb']:
    makefile_options['PREPROCESSOR_FLAGS'] += '-DABSORB '

makefile_options['PREPROCESSOR_FLAGS'] += '-DNGHOST=' + str(args['nghosts']) + ' '

# Step 3. Create new files, finish up
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)
makefile_template = re.sub('# Template for ', '# ', makefile_template)
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

# Finish with diagnostic output
print('==============================================================================')
print('Your TRISTAN-FF distribution has now been configured with the following options:')
if (specific_cluster):
    if (args['perseus']):
        print('  Cluster configurations:  `Perseus`' )

print('SETUP ........................................................................')
print('  Userfile:                ' + makefile_options['USER_FILE'])
print('  # of ghost zones:        ' + str(args['nghosts']))
print('  Absorbing boundaries:    ' + ('ON' if args['absorb'] else 'OFF'))

# print('PHYSICS ......................................................................')
# print('  External fields:         ' + ('ON' if args['extfields'] else 'OFF'))
# print('  Absorbing boundaries:    ' + ('ON' if args['absorb'] else 'OFF'))
# print('  Cooling:                 ' + args['radiation'])
# print('  Photon emission          ' + ('ON' if args['emit'] else 'OFF'))
# print('  QED step                 ' + ('ON' if args['qed'] else 'OFF'))
# print('  BW pair production       ' + ('ON' if args['bwpp'] else 'OFF'))

print('TECHNICAL ....................................................................')
print('  Compiler:                ' + ('intel' if args['intel'] else 'gcc'))
print('  Debug mode:              ' + ('ON' if args['debug'] else 'OFF'))
print('  Output:                  ' + ('HDF5' if args['hdf5'] else 'binary'))
print('  MPI version:             ' + ('old' if not args['mpi08'] else 'MPI_08'))
print('  `IFPORT` mkdir:          ' + ('ON' if args['ifport'] else 'OFF'))

print('==============================================================================')

print('  Compilation command:     ' + makefile_options['COMPILER_COMMAND'] \
    + makefile_options['PREPROCESSOR_FLAGS'] + makefile_options['COMPILER_FLAGS'])
