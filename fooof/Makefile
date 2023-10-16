.PHONY: ALWAYS
id_age_sex_ses.csv: ALWAYS
	./00_id_age_sex_ses.R | mkifdiff $@

python/fooof/Results/allSubjectsFooofMeasures_20230516.csv: $(wildcard /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Resting_State/AfterWhole/ICAwholeClean_homogenize/1*_2*_*_Rem_rerefwhole_ICA_icapru.set)
	python python/fooof/runFooof.py
