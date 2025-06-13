I had issues with the DRAM (the first DRAM, not DRAM2) database set up, and the first thing that some people say it will fix it is by modifying the script:
/home/daniel.castanedamogo/anaconda3/envs/DRAM1.2/lib/python3.11/site-packages/mag_annotator/database_processing.py

In here I had to change this:

    `merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)`

To this


    `merge_files(glob(path.join(hmm_dir, 'hmm', 'VOG*.hmm')), vog_hmms)`
