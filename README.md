# n-body-gravity
A simulation of the n body problem involving bodies of equal mass orbiting under gravity

This is a scalable simulation of gravity. I am using rescaled units so if you want to use this for real data you will need to find a conversion to "real" unit. 

If you want to scale this simulator you can go to the ```#define N 3``` and change it to the number of dimentions you want to simulate. 

You can also change ```#define D 3``` to change the number of dimentions (currently it is in 3d). I wouldn't reccomend going to one dimentions, and higher than 3 would be hard to visualize, but you do you. 

My plotting sceme uses Mathematica which is a propriatary program this will work with any number of bodies. There is also a python program that will animate for three bodies in three dimentions specifically. I am working to generalize this program to n bodies. The python animation also shows the path of the bodies as the animation progresses.  
