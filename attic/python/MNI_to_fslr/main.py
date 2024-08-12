# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

from neuromaps.datasets import fetch_annotation
from neuromaps import transforms
neurosynth = fetch_annotation(source='neurosynth')
fslr = transforms.mni152_to_fslr(electrodeLocs_MNIcoordinates_cortex_atlas.nii.gz, '32k')