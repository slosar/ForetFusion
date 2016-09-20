mode: bycore
N: 4
threads: 1
hostfile: auto
job_name: sdss_catalog
command: |
     source ~/.bashrc;
     OMP_NUM_THREADS=%threads% mpirun -hostfile %hostfile% python run_catalog.py >test_catalog.log 2>test_catalog.err
