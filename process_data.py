import read_GCAM_DB
import process_GCAM_data
import plotting_script

if __name__ == '__main__':
    ### double check constants.py before running this code to ensure correct file placement and contents ###
    read_GCAM_DB.main()
    process_GCAM_data.main()
    plotting_script.main()
