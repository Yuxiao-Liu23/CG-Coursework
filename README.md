Hi，here is the README.md for CG-Coursework (Finally !!!)

So its a renderer written in C++17,
1. What function it have?
Here is some function i made:

OBJ geometry and material file loading 
Wireframe 3D scene rendering 
Flat colour 3D scene rasterising 
Keyboard control of camera position 

Frame saving and video compositing (I will explain it later)
Scene contains more than one model (Cornell box, ball for goroud/phong,and logo on /models and /resources)(no eiffel,its just a test if you want to add more model in the future...... so ignore it)
Use of depth buffer to resolve occlusion
Keyboard control of camera orientation //

Hard Shadow (without soft edges) plus ambient lighting
Some form of surface texture mapping 
Diffuse lighting 
Simple animation (orbit, lookAt, fly-through,......)

Gouraud shading
Rough attempt at soft shadows 
Mirrored surfaces
Specular lighting 

Phong Shading
Smooth and elegant soft shadows(maybe lol)

2.Some thing to say?

I'd like to explain the implemented functionality to help you operate this renderer:
First, I set a Boolean variable `Mode` for each mode and toggled them using buttons (e.g., wireframe mode...). I think this might make editing easier? If you want to implement camera rotation, disabling Lookat mode might work better (because the lens won't always be tracking the object) (I remember lookat is enabled by default lol).

Oh, I also want to talk about soft shadows. If you want soft shadows to be more obvious, you can increase the number of light sources (8). However, be careful not to make the number too large, otherwise the program may crash!

3.Button binding!!! (which is basically the instruction manual: to be honest, without this I would forget the function of some buttons, haha):

0	Ray Tracing Mode
1	Wireframe Rendering
2	Flat Rasterised Rendering
3	Proximity Lighting
4	Diffuse Lighting
5	Gouraud Shading
E	Phong Shading (per-pixel)
6	Toggle Texture Mapping

W / S	Move camera forward / backward (Z axis)
← / →	Move camera left / right (X axis)
↑ / ↓	Move camera up / down (Y axis)

A	Orbit left (Y axis +)
D	Orbit right (Y axis –)
R	Orbit upward (X axis +)
F	Orbit downward (X axis –)

J	Rotate yaw + (Y axis) 
L	Rotate yaw – 
I	Rotate pitch + (X axis) 
K	Rotate pitch – 

O	Toggle automatic orbit
P	Reset view using lookAt(sceneCentre)

T / G		Move light up / down (Y axis)
Y / H		Move light right / left (X axis)
U / M	Move light forward / backward (Z axis)
9		Toggle Hard Shadow ↔ Soft Shadow

7	Load Cornell Box scene
8	Load Sphere model
Z	Insert Logo model into Cornell Box

V	Toggle frame recording ON/OFF

IMPORTANT !!!!!:

gouraud and phong might perform better in ball model.but for phong it fine because i use Specular light for Cornell model and Phong light for Ball

For key v I'd like to add something. I designed this feature to automatically capture each frame and output a PPM file. Initially, it worked in conjunction with a mode called Demo mode. Demo mode was designed to automatically plan what the renderer should do at each stage (e.g., 0-200 frames: wireframe rendering with automatic rotation; 201-400 frames: rasterization without rotation; 401-600 frames: switch to a different model, etc.). Basically, you could customize any mode you wanted and let them run automatically, as long as you defined the frame range. However, I later discovered significant limitations: for example, when automatically running to 200 frames, although it might switch from rasterization to ray tracing, the background was completely obscured, making it very difficult to control the frame count. Manual operation was much better for controlling the details (you know when to switch for the best effect!). So I set it to false (not enabled). However, I still kept all the demo mode code; having more options is always good, haha.

4.How to run it?
Okay, I admit this is crucial:

If everything is working correctly, you need to run `make speedy` to display the correct speed (because the computation is really huge!).

Additionally, I designed `make Cornell-movie` to handle frame saving and video output, but I'm still debugging it so I'm not sure lol.

Since my code was written and compiled on Windows, it sometimes throws errors on Linux (for example, "infinity" works perfectly on Windows, but crashes and fails to display images on Linux—probably due to a parameter in speedy—I eventually solved the problem by replacing infinity with max()).

5. In general, due to some differences between Windows and Linux, there may be some unknown errors (I tried testing on an MVB computer, but why do so many computers have problems preventing the test?). I encountered an Infinity error when I was fortunate enough to successfully turn on my computer before (haha, just kidding).
