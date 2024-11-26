
# Ray Tracing

In this [assignment](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/87d3b4c5ba0c41f34991e932270fcaf08d0b4cf6/Offline-3-ray-casting/CSE%20410%20Offline%203%20Specs.pdf), I have with OpenGL library to generate realistic images for a few geometric shapes using ray tracing with appropriate illumination techniques.



## Features

- **Controlling Camera**
  - Translation
    - up arrow –move forward
    - down arrow –move backward
    - right arrow –move right
    - left arrow –move left
    - page up –move up
    - page down –move down
  - Rotation
    - 1 –rotate/look left
    - 2 –rotate/look right
    - 3 –look up
    - 4 –look down
    - 5 –tilt counterclockwise
    - 6 –tilt clockwise
- Called a function named ```loadData()``` from my main function. This function reads a text file named ```scene.txt``` (input file) containing details of different objects (shapes) and lights present in the environment. 
- There can be three types of objects (shapes) in the file - sphere, triangle, and object (shape) having a general quadratic equation in 3D.
- There will be two types of light sources - point lights and spotlights. In the case of point lights, their color and position will be specified. In the case of the spotlights, their direction and cutoff angle will be specified in addition to the position and color information. 
- The complete spec is [here](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/87d3b4c5ba0c41f34991e932270fcaf08d0b4cf6/Offline-3-ray-casting/CSE%20410%20Offline%203%20Specs.pdf).
- Created a separate header file named ```1905003_Header.h```. Included this file in the cpp file containing the main function. 
- Each time 0 is pressed in the *OpenGL* window (NOT Codelbocks console), a ```.bmp``` image is generated and stored in the same folder as ```main.cpp```, e.g., ```Output_11.bmp```, ```Output_12.bmp``` and so on.




## Language, Libraries and IDE

**Programming Language:** C++

**Library:** OpenGL

**IDE:** Codeblocks


## Run Locally

- Follow the instructions for setting up OpenGL in Codeblocks for Windows in [Setting Up OpenGL.docx](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/80c3b76d5e1fe36637a03e311270df7e70f4e223/Offline-1-OpenGL/Setting%20Up%20OpenGL.docx)
- Upon creating a new project in Codeblocks, a ```main.cpp``` file will be automatically created.
- Copy the following the files from [this](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/tree/master/Offline-3-ray-casting) folder into the folder where the project is located:
  - ```1905003_Header.h```
  - ```bitmap_image.cpp```
  - ```bitmap_image.hpp```
  - ```scene.txt``` (can be replaced with contents from ```test-input.txt``` for creating a different scene)
- Copy-paste contents of ```1905003_Main.cpp``` from the folder into ```main.cpp```
- ```camera-locations.txt``` contains suggestions for good camera locations (optional to use in ```main.cpp```)
- Click *Build and run* button at the top of the IDE
- The output files ```Output_11.bmp```, ```Output_12.bmp``` etc image will be created in the same folder as ```main.cpp```.


## Demo

https://github.com/user-attachments/assets/6f337853-9556-4cc2-b87e-107b09f25ed2

## Author

- [Sheikh Intiser Uddin Dipta](https://github.com/sheikhDipta003)

