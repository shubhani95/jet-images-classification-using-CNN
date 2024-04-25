# jet-images-classification-using-CNN

In this study, we seek to utilize ML techniques to build a classifier for finding signs of new physics at the LHC using 2HDM Type-II. In particular, we will incorporate the advanced technique of image recognition by designing a CNN and visualizing what is ‘seen’ in a detector in particle physics experiments. Of course, there are well-established LHC searches using traditional techniques and alternative clustering algorithms for evaluating 2HDM final states. Here, we will deploy an image recognition-informed jet tagger’s exceptional ability to map the jet-level information to an image and distinguish signals from relevant backgrounds.

We will employ jet-level image recognition studies. In the context of boosted decays such that the entire signal event can be clustered into a single fat jet using a large cone size, and the big discriminatory feature for signal and background is the presence of a two-prong jet substructure, as done in the context of HSM → bb~ decay. For our study, we focus on b-jet final states from 2HDM Type-II with decay chains of the form gg →H→hh→bb~bb~. We also consider leading SM backgrounds, such as: gg,qq'→bb~bb~, gg, qq' → Zbb~, and gg, qq' → tt~.

# thesis_sjain_images_classification_notes.pdf
This document encompasses project notes, detailing tools, event generation, analysis cutflow, results, and future prospects.

# event_generation.txt
Contained within this file are comprehensive details regarding event generation and parton-level cuts to be implemented in Madgraph, encapsulating essential guidelines for precise simulation and analysis.

# ML_fatjetimages_delphes.cpp
Enclosed herein is the expert mode C++ code of Madanalysis5, meticulously tailored for analysis purposes, facilitating the extraction of essential information and variables required for generating .csv files. These files serve as crucial inputs for crafting jet images indispensable for classification tasks.

# pysaf.py
This Python file serves as a versatile tool for plotting the distributions extracted from the final .SAF output of Madanalysis5.

# pre_processing.ipynb
This Jupyter notebook is crafted for preprocessing .csv data, crafting jet images, and orchestrating the entire pipeline of training and deploying a convolutional neural network (CNN) model for image classification. It seamlessly integrates data preprocessing, model training, evaluation, and performance analysis, offering insights into CNN outputs and performance metrics crucial for classification tasks.

# distribution_plots_fixedr_vs_varr
This folder contains kinematic distributions obtained through the Anti-kt algorithm, employing a fixed-R value of 1.2 and variable-r implementation. These distributions offer valuable insights into the clustering behavior, facilitating a comprehensive understanding of the underlying physics phenomena.

# avg_jet_images
Contained within this folder are representative average jet images for both signal and background events. These images offer a succinct portrayal of the collective features of N events, mitigating the impracticality of navigating through a vast gallery. By presenting these averaged representations, we provide a concise yet informative snapshot of the jet characteristics crucial for discerning the signal from backgrounds.

# cnn_results
This folder houses critical performance metrics, including accuracy plots, AUC curves, and CNN scores, providing comprehensive insights into the classification efficacy and model performance.

# npz_h5_files
This folder contains essential .npz files indispensable for the image classification task and training of the CNN model. 










