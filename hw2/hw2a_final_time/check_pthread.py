import os
import sys
import json
import hashlib
import subprocess

# Check Fast
for i in range(1, 17):
  testcase = str(i)
  if i < 10:
    testcase = '0' + testcase
  if len(sys.argv) > 1:
    testcase = str(sys.argv[1])
  testcase = "fast" + testcase


  with open(f'../testcases/{testcase}.txt', 'r') as f:
    line = f.readline()
    parm_list = line.split(" ")
    print('\033[1;34;48m' + '='*35 + f' Start running test {testcase} ' + '='*35 + '\033[1;37;0m')

    out_file_name = f'./out.png'
    os.system(f'rm {out_file_name}')

    main_script = f'srun -n1 -c4 ./hw2a {out_file_name} {parm_list[0]} {parm_list[1]} {parm_list[2]} {parm_list[3]} {parm_list[4]} {parm_list[5]} {parm_list[6]}'
    print(main_script)
    os.system(main_script)
    print('done.')

    print('')
    print('verifying ...')
    verify_script = f'hw2-diff ../testcases/{testcase}.png {out_file_name}'
    print(verify_script)
    os.system(verify_script)

# ------------------------------------------------------------------------------------------------------------------------ 

# Check Slow 
for i in range(1, 17):
  testcase = str(i)
  if i < 10:
    testcase = '0' + testcase
  if len(sys.argv) > 1:
    testcase = str(sys.argv[1])
  testcase = "slow" + testcase


  with open(f'../testcases/{testcase}.txt', 'r') as f:
    line = f.readline()
    parm_list = line.split(" ")
    print('\033[1;34;48m' + '='*35 + f' Start running test {testcase} ' + '='*35 + '\033[1;37;0m')

    out_file_name = f'./out.png'
    os.system(f'rm {out_file_name}')

    main_script = f'srun -n1 -c4 ./hw2a {out_file_name} {parm_list[0]} {parm_list[1]} {parm_list[2]} {parm_list[3]} {parm_list[4]} {parm_list[5]} {parm_list[6]}'
    print(main_script)
    os.system(main_script)
    print('done.')

    print('')
    print('verifying ...')
    verify_script = f'hw2-diff ../testcases/{testcase}.png {out_file_name}'
    print(verify_script)
    os.system(verify_script)


# ------------------------------------------------------------------------------------------------------------------------ 

# Check Strict 
for i in range(1, 37):
  testcase = str(i)
  if i < 10:
    testcase = '0' + testcase
  if len(sys.argv) > 1:
    testcase = str(sys.argv[1])
  testcase = "strict" + testcase


  with open(f'../testcases/{testcase}.txt', 'r') as f:
    line = f.readline()
    parm_list = line.split(" ")
    print('\033[1;34;48m' + '='*35 + f' Start running test {testcase} ' + '='*35 + '\033[1;37;0m')

    out_file_name = f'./out.png'
    os.system(f'rm {out_file_name}')

    main_script = f'srun -n1 -c4 ./hw2a {out_file_name} {parm_list[0]} {parm_list[1]} {parm_list[2]} {parm_list[3]} {parm_list[4]} {parm_list[5]} {parm_list[6]}'
    print(main_script)
    os.system(main_script)
    print('done.')

    print('')
    print('verifying ...')
    verify_script = f'hw2-diff ../testcases/{testcase}.png {out_file_name}'
    print(verify_script)
    os.system(verify_script)