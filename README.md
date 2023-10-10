# replication_bstar
Replication of the bSTAR sequence and open-source implementation (Pulseq + BART)

## General use

Install BART in your system and add modified lines from this repo to the original BART source code: https://github.com/namgyunlee/bart_para_fista

For Windows users, BART must be installed in Windows WSL.

For Pulseq bSTAR acquired on Siemens scanners, convert TWIX format to ISMRMRD format using this m-file: `demo_convert_siemens_to_ismrmrd.m`.

Make a new 'json' file for each dataset. Here is an example of a 'json' file:

	{
       "siemens_twix_file": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/raw/meas_MID00539_FID62952_pulseq_bstar_1_38ms_1_61mm_bh_i4_17k.dat",
       "ismrmrd_data_file": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/raw/meas_MID00539_FID62952_pulseq_bstar_1_38ms_1_61mm_bh_i4_17k.h5",
      "ismrmrd_noise_file": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/raw/noise_meas_MID00539_FID62952_pulseq_bstar_1_38ms_1_61mm_bh_i4_17k.h5",
	            "seq_file": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/seq/bstar_ecg0_TR1.38ms_1.61mm_b1929_rf200_i4_17k_FA25_self0_WASP.seq",
                "trj_file": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/trajectory/traj_b1929_n288_Zhao_2020_MRM_v2_fine_tuning.trj",
             "output_path": "D:/data_pulseq_bstar/pulseq_bSTAR_USC_ACR_phantom_20230406/pulseq_bstar_1_38ms_1_61mm_bh_i4_17k",
               "bart_path": "/home/image/bart",
              "study_date": "20230406",
      "recon_parameters": [
        {
	      "recon_matrix_size": [ 360, 360, 360 ],
          "dicom_matrix_size": [ 320, 320, 320 ],
          "recon_interp_factor": 1.060465416011151,
          "cal_size": [ 32, 32, 32 ],
          "lambda": 0.0001,
          "max_iter": 30,
          "B": 1,
          "kappa": 1
        }
      ],
      "noir_conf": [
        {
	      "a": 16,
          "b": 16,
          "max_iter": 25
        }
      ]
	}


Update variables in `batch_pulseq_bstar_acr_phantom_20230614.m`.

    package_directory = 'path-to-this-package';
    ismrmrd_directory = 'path-to-ISMRMRD-package';
    pulseq_directory  = 'path-to-pulseq-package';

Run `batch_pulseq_bstar_acr_phantom_20230614.m`.

## Free-breathing Pulseq bSTAR (code not available yet)

In `demo_estimate_resp.m`, I used source code from 
1. Dr.Frank Ong's extreme_mri Python package: https://github.com/mikgroup/extreme_mri
2. Dr.Li Feng's XD-GRASP MATLAB package: https://cai2r.net/resources/xd-grasp-matlab-code