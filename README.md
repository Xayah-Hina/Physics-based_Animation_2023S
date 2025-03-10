# Physics-based Animation 4860-1081 2023S

## Task Results (native implemetation)

|             task 1             |             task 2             |             task 3             |             task 4             |              task 5              |
|:------------------------------:|:------------------------------:|:------------------------------:|:------------------------------:|:--------------------------------:|
| ![img_1.png](images/img_1.png) | ![img_2.png](images/img_2.png) | ![img_3.png](images/img_3.png) | ![img_4.png](images/img_4.png) |  ![img_5.png](images/img_5.png)  |
|             task 6             |             task 7             |             task 8             |             task 9             |             task 10              |
| ![img_6.png](images/img_6.png) | ![img_7.png](images/img_7.png) | ![img_8.png](images/img_8.png) | ![img_9.png](images/img_9.png) | ![img_10.png](images/img_10.png) |

## Task Results (houdini implemetation)

For better visualization results and parameters tuning, I implemented the tasks in Houdini.

### Build Instructions

To build Houdini HDK Node, you need to have Houdini installed on your system. You can download the free version of Houdini from [here](https://www.sidefx.com/products/houdini-apprentice/).
After houdini is installed, you need to set the `Houdini_PATH` in the `CMakeLists.txt` file to the path where Houdini is installed.
```cmake
set(Houdini_PATH "C:/Program Files/Side Effects Software/Houdini 20.5.487") # replace with your Houdini installation path
```
then run the following commands to build the project.
```bash
cmake -B build -S .
cmake --build build
```
After building the project, you can open the project files and find the compiled HDK node.