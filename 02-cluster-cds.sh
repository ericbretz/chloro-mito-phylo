#!/usr/bin/env bash

rm -rf ./cp/*
mkdir -p ./cp
vsearch --cluster_fast ./cp.fasta --clusterout_id --clusters ./cp/cp_cluster_ --centroids ./cp_centroids.fasta --threads 14 --iddef 1 --id 0.5 &> ./cp.log

rm -rf ./mt/*
mkdir -p ./mt
vsearch --cluster_fast ./mt.fasta --clusterout_id --clusters ./mt/mt_cluster_ --centroids ./mt_centroids.fasta --threads 14 --iddef 1 --id 0.5 &> ./mt.log
