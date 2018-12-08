# MoG-for-small-target-detection
This source code is for the infrared small target detection method based on  Mixture of Gaussians (MoG) with Markov random field (MRF) proposed in our paper of PR2018: Chenqiang Gao, Lan Wang, Yongxing Xiao, Qian Zhao,Deyu Meng, “Infrared small-dim target detection based on Markov random field guided noise modeling, Pattern Recognition, vol. 76, no. Supplement C, pp. 463-475, 2018/04/01/, 2018.  
Note that this source code does not contain the threhold segmentation module, which can be easily implemented.

If you use this code in your publications, please cite:  
@article{Gao2018Infrared,  
   author = {Gao, Chenqiang and Wang, Lan and Xiao, Yongxing and Zhao, Qian and Meng, Deyu},  
   title = {Infrared small-dim target detection based on Markov random field guided noise modeling},  
   journal = {Pattern Recognition},  
   volume = {76},  
   number = {Supplement C},  
   pages = {463-475},  
   month = {2018/04/01/},  
   year = {2018}  
}  

If you use the test images in your publications, in additoin to above reference, please cite the following references, too:

@article{Gao2013IPI,  
   author = {Gao, Chenqiang and Meng, Deyu and Yang, Yi and Wang, Yongtao and Zhou, Xiaofang and Hauptmann, Alex},  
   title = {Infrared Patch-Image Model for Small Target Detection in a Single Image},  
   journal = {Image Processing, IEEE Transactions on},  
   volume = {22},  
   number = {12},  
   pages = {4996-5009},  
   year = {2013}  
}

@article{Gao2012,  
 author = {Chenqiang, Gao and Tianqi, Zhang and Qiang, Li},  
 title = {Small infrared target detection using sparse ring representation},  
  journal = {IEEE Aerospace and Electronic Systems Magazine},  
   volume = {27},  
   number = {3},  
   pages = {21-30},  
   year = {2012}  
}  


## How to use this code?
You just run the the file of *main.m* in matlab. Our matlab version for development is Matlab2018a. 
Some modules are written by C++, and you should recompile them if you run
it with errors, Like as follows:
> mex update_Z.cpp  
> mex reconstructImage.cpp  

## Contact
If you have any questions, please contact:  
Author: Chenqiang Gao  
Email: gaochenqiang@gmail.com *or* gaocq@cqupt.edu.cn  
Copyright: Chongqing University of Posts and Telecommunications  
## License
This code is only freely available for non-commercial research use, but please contact us if for other purposes.

If you find some help for you, star is a good reward ^_^. 