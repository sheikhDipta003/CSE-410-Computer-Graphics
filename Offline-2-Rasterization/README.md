
# Raster-based graphics pipeline

In this [assignment](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/370745e19e9a7285f90727e9b1141cc3c5ee3825/Offline-2-Rasterization/Offline%202%20Specifications.pdf), I have developerd the raster-based graphics pipeline used in OpenGL. The pipeline can be thought of as a series of six stages. I implemented roughly 4 stages of the pipeline.
- Stage 1: modeling transformation
- Stage 2: view transformation
- Stage 3: projection transformation
- Stage 4: clipping & scan conversion using Z-buffer algorithm

## Features

- There are two input files: *scene.txt* and *config.txt*. The first file will contain the scene description and second file will contain the necessary information for Z-buffer algorithm.
- *scene.txt* contains the following lines:
  - Line 1: eyeX eyeY eyeZ
  - Line 2: lookX lookY lookZ
  - Line 3: upX upY upZ
  - Line 4: fovY aspectRatio near far
  - Lines 1-3 of scene.txt state the parameters of the gluLookAt function and Line 4 provides the gluPerspective parameters.
  - The subsequent lines contain a command followed by the parameters necessary to execute that command in OpenGL. See details in the [spec](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/370745e19e9a7285f90727e9b1141cc3c5ee3825/Offline-2-Rasterization/Offline%202%20Specifications.pdf).
- *config.txt*'s  first Line contains two integers, representing ```Screen_Width``` and ```Screen_Height``` respectively. The final rendered image will have the width and height equal to these values respectively.
- This program outputs five files: *stage1.txt, stage2.txt, stage3.txt, z-buffer.txt, and out.bmp*. The first three files will contain the output of the first three stages, respectively. The fourth file will contain z-buffer values (only those which are less than the max value). And the fifth file will be a bmp image generated by the pipeline.




## Language, Libraries and IDE

**Programming Language:** C++

**Library:** OpenGL

**IDE:** Codeblocks


## Run Locally

- Follow the instructions for setting up OpenGL in Codeblocks for Windows in [Setting Up OpenGL.docx](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/80c3b76d5e1fe36637a03e311270df7e70f4e223/Offline-1-OpenGL/Setting%20Up%20OpenGL.docx)
- Upon creating a new project in Codeblocks, a ```main.cpp``` file will be automatically created.
- Copy all the files from [this](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/tree/370745e19e9a7285f90727e9b1141cc3c5ee3825/Offline-2-Rasterization) folder into the folder where the project is located. You may need to manually copy-paste contents of ```main.cpp``` from the folder into the project's default ```main.cpp```
- In the first line of the ```main()``` function, insert the filepath for ```scene.txt``` from any folder (```1``` ,``` 2``` ,``` 3``` ,``` 4``` ,``` 5```) under *IOs* folder 
- Click *Build and run* button at the top of the IDE
- The output files and ```out.bmp``` image will be created in the same folder as ```main.cpp```.


## Demo


https://github.com/user-attachments/assets/f54547dd-14fe-4406-be3f-26514fa7f77a


## Author

- [Sheikh Intiser Uddin Dipta](https://github.com/sheikhDipta003)

