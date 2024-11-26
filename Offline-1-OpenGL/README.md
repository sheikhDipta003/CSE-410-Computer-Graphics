
# Rolling Ball and Magic Cube with Controllable Camera

Code for assignment 1 of CSE410 sessional course. In this assignment, a keyboard-contraollable camera is implemented, along with a rolling ball on the xy plane and a magic cube that can be transformed back and forth between an octahedron and a sphere.

## Features

- **Task 1 : Controlling Camera**
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
- **Task 2: Rolling Ball**
  - A ball is drawn on the xy plane. It has a direction (in the xy plane) and a movement speed.
  - The ball is drawn by dividing it into sectors and stacks. Neighbouring sectors have alternating colours.
  - An arrow is shown in the current direction of the ball.
  - The floor is drawn as a checkerboard, change the length and width of the checkerboard in the code.
  - There are two modes: manual control and simulation. To shift between modes, press *comma ( , )*.
    - *Manual Control*
      - j –direction will rotate counterclockwise
      - l –direction will rotate clockwise
      - i –go forward
      - k –go backward
    - *Simulation mode*: In this mode, I implemented time-driven simulation by considering time in discrete units. The ball moves in a region bounded by planes (not necessary 4). The normals of the planes are coplanar with the balls' direction. When a collision occurs, the direction of the ball is updated using the reflection formula. 
      - j –direction will rotate counterclockwise
      - l –direction will rotate clockwise
- **Task 3 : Magic Cube**
  - The magic cube makes a transition between a sphere and an octahedron.
  - Pressing the following keys will change the shape:
    - . (dot) – sphere to octahedron
    - , (comma) – octahedron to sphere
  - Only one triangle, one cylinder segment and one sphere segment are drawn, and then various transformations - translation, rotation, scaling are used to put them in the right places.



## Language, Libraries and IDE

**Programming Language:** C++

**Library:** OpenGL

**IDE:** Codeblocks


## Run Locally

- Follow the instructions for setting up OpenGL in Codeblocks for Windows in [Setting Up OpenGL.docx](https://github.com/sheikhDipta003/CSE-410-Computer-Graphics/blob/80c3b76d5e1fe36637a03e311270df7e70f4e223/Offline-1-OpenGL/Setting%20Up%20OpenGL.docx)
- Upon creating a new project in Codeblocks, a ```main.cpp``` file will be automatically created.
- Task-2:
  - Copy contents of the file ```rolling_ball.cpp``` into ```main.cpp```, click *Build and run* button at the top of the IDE
- Task-3:
  - Copy contents of the file ```magic_cube.cpp``` into ```main.cpp```, click *Build and run* button at the top of the IDE

## Demo

Insert gif or link to demo


## Author

- [Sheikh Intiser Uddin Dipta](https://github.com/sheikhDipta003)

