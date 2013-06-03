Generate list of tests
======================

Generate list of commands for data generation and sampling,
for values of sigma0, ntrials, and ndata specified in `generate_command.py`,
and for all complexes listed in `complex_list_3` and `complex_list_4`.

## Usage:

`python generate_command.py complex_list_3 complex_list_4`

The files `GENERATE_COMMAND_LIST` and `SAMPLE_COMMAND_LIST` can be created with

`python generate_command.py complex_list_3 complex_list_4 | grep generate > GENERATE_COMMAND_LIST`
`python generate_command.py complex_list_3 complex_list_4 | grep sample   > SAMPLE_COMMAND_LIST` 
