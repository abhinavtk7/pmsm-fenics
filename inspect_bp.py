'''
Python script to inspect .bp file and its variables
'''

import adios2

# Open the .bp file
with adios2.open("/root/shared/pm3Dec30_TEAM30_three_50.0_1_100_False/Az.bp", "r") as fh:
    for step in fh:
        print("Step:", step.current_step())
        for name in step.available_variables():
            print("Variable:", name)
            print("Metadata:", step.available_variables()[name])
