# DNA Sequence Analysis Using the MinHash Algorithm

This repository contains the code and documentation for a project focused on the efficient comparison and analysis of large sets of DNA sequences using the MinHash algorithm. The primary goal of this project is to align short DNA reads to a reference genome and identify the most similar regions (bins) within the reference.

## Overview

In this project, we process and analyze genomic data from various organisms to gain insights into their genetic makeup, evolutionary history, and functional characteristics. We use the MinHash algorithm to efficiently compare DNA sequences and reduce the feature dimension, enabling faster and more effective analysis of large genomic datasets.

## Usage

The main script for this project is Milestone_Script.py. The script will read the input FastQ and Fasta files containing DNA reads and reference genomes, respectively. It will then perform the MinHash algorithm to generate MinHash signatures for each DNA read and reference bin, followed by calculating the Jaccard similarity between them. The results will be stored in a list, which can be used for further analysis or visualization.

## Results and Evaluation

The MinHash algorithm effectively reduces the dimensionality of the DNA sequence data, enabling efficient comparisons between the reads and reference bins. The results show a strong correlation with the benchmark dataset, indicating that our approach successfully aligns short DNA reads to the reference genome and identifies the most similar regions. This method can be further refined and optimized for various applications in genomics, such as identifying structural variants or studying the evolution of genomes.

## Contributors

Yingxiang (Andy) Ma

Peter Tian
