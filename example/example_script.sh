
cd 1d7p

../helper_scripts/gromacs_psi_phi_chi_max 1d7p.prmtop

cd gromacs_psi_phi_chi_max

# rename files because numbering of residues differs from 1d7p to 1iqd
../../helper_scripts/rename_filename_residues.py

cd ../..

cd 1iqd_C2

../helper_scripts/gromacs_psi_phi_chi_max 1iqd_C2.prmtop

cd ..

../KBL.py 1d7p/gromacs_psi_phi_chi_max 1iqd_C2/gromacs_psi_phi_chi_max -a phi,psi -r 20,21,22,23,24,25,26,27,28,29,30 -suff "hairpin"


# now, with pymol one should open 1iqd_C2_WT/1iqd.pdb and load the file with the KBL values kbl_1d7p_1iqd_C2_phi_psi_hairpin.pymol

# we also can further investigate the behaviour of backbone angles throughout the trajectory with time_scale_hist.py

../time_scale_hist.py -f 1iqd_C2/gromacs_psi_phi_chi_max/psiPHE27.xvg -pre "1iqd_C2 " # creates file psiPHE27.png in folder 1iqd_C2/gromacs_psi_phi_chi_max/

../time_scale_hist.py -f 1d7p/gromacs_psi_phi_chi_max/psiPHE27.xvg -pre "1d7p " # creates file psiPHE27.png in folder 1d7p/gromacs_psi_phi_chi_max/

