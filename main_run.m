%   ==========
%   The Code is created based on the method described in the following paper 
%   "Computed Tomography Ring Artifact Correction Method with Super-Pixel and Adaptive Relative Total Variation", Na Li.
%   The paper was under review in the journal
%   The code and the algorithm are for non-comercial use only.
%  
%   Na Li (na.li3@gdmu.edu.cn)
%   Date  : 17 Feb 2024
%   Version : 1.0 
%   Copyright 2024, Guangdong Medical University.

clear
close all
clc

win_width=1600;
win_lev=400;

load('ring_img_01.mat');

[no_ring_img_t] = ring_remove(img);
artifact = double(img)-double(no_ring_img_t);
figure,imshow(img(:,:,round(size(img,3)./2)),[win_lev-(win_width./2) win_lev+(win_width./2)]);
figure,imshow(no_ring_img_t(:,:,round(size(img,3)./2)),[win_lev-(win_width./2) win_lev+(win_width./2)]);
figure,imshow(artifact(:,:,round(size(img,3)./2)),[-250 250]);
