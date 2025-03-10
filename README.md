# YOLOv8 Autofocus for Droplet Imaging

## Overview

This project implements an automated droplet imaging system using **YOLOv8-based autofocus** and **Hough Transform-based drift correction**. The script dynamically finds the focal plane and corrects for **Z-drift** and **XY-drift**, ensuring high-precision imaging over time.

## Features

- **Deep Learning-Based Autofocus**: Uses a **YOLOv8 model** to determine whether an image is in focus and adjust the Z-position accordingly.
- **XY-Drift Correction**: Uses **Hough Transform** to detect droplet centers and correct movement.
- **Time-Lapse Imaging**: Captures a sequence of images at defined time intervals.
- **Z-Stack Image Acquisition**: Collects images at multiple focal depths for each droplet.
- **Logging & Tracking**: Saves focus history in a CSV file for analysis.

## Requirements

### **Software Dependencies**

Ensure the following Python libraries are installed:

```bash
pip install numpy pandas opencv-python tisgrabber ultralytics serial
```

### **Hardware Requirements**

- Microscope with **motorized stage**
- Camera compatible with **TIS Imaging Control (IC)**
- **PriorStage motorized stage controller**

## Installation

1. **Clone the repository**:

```bash
git clone https://github.com/GielenLab/YOLOv5-v8-polymicrobial/YOLOv8-Autofocus.git
cd YOLOv8-Autofocus
```

2. **Ensure hardware is connected** (camera & motorized stage).

## Usage

### **Running the Script**

Execute the Python script:

```bash
python xyz_focus_YOLOv8.py
```

### **Configurable Parameters**

Modify the parameters in the script for customized imaging:

| Parameter       | Description                        | Default Value |
| --------------- | ---------------------------------- | ------------- |
| `Width`         | Image width                        | 640           |
| `Height`        | Image height                       | 480           |
| `total_nb`      | Number of time points              | 500           |
| `time_interval` | Time between images (sec)          | 45            |
| `focus_step`    | Step size for focus (0.1um)        | 8             |
| `nb_slices`     | Number of slices above/below focus | 10            |

### **Focus History Logging**

The script records autofocus behavior in `focus_history.csv`, storing:

- Timepoint
- Timestamp
- Droplet ID
- Focus Level
- Confidence Score
- Adjustments made

## How It Works

### **Autofocus Process**

1. **YOLOv8 model** classifies focus levels (Above, Below, InFocus).
2. **Stepwise adjustments** are made to reach the focal plane.
3. **If misalignment is detected**, it dynamically corrects Z-position.

### **XY-Drift Correction**

1. **Hough Transform** detects droplet position.
2. If drift is found, the **motorized stage** moves to recenter the droplet.

### **Image Acquisition**

- Multiple slices are taken **above and below the focal plane** to create a Z-stack.
- Images are saved in **separate folders** for each droplet.

# **MATLAB Scripts for Image Analysis**

The repository also includes MATLAB scripts for bacterial image processing:

## `extract_clusters_dual_model_date_from_txt_files.m`

- Extracts and processes bacterial cluster information from detection files.
- Identifies **in-plane and out-of-plane cells** and corrects for motion.
- Saves the results as a video and an Excel file for further analysis.

## `extract_clusters_from_text_files.m`

- A streamlined version of the cluster extraction script.
- Detects bacterial clusters and logs cell counts over time.
- Saves the processed data to an Excel file.

## `video_with_graph_1class.m`

- Processes images and overlays bacterial count graphs in real time.
- Displays both **timelapse images** and **corresponding cell count plots**.
- Saves the output as a video.

## `video_with_graph_2class.m`

- Similar to `video_with_graph_1class.m`, but distinguishes between **two bacterial classes** (e.g., **Pseudomonas** and **Staphylococcus**).
- Saves annotated videos with bacterial count graphs.

# Troubleshooting

## **Common Issues & Solutions**

| Issue                      | Possible Cause               | Solution                            |
| -------------------------- | ---------------------------- | ----------------------------------- |
| No images captured         | Camera not initialized       | Ensure camera is properly connected |
| Focus not adjusting        | YOLO model misclassification | Train or fine-tune the model        |
| Incorrect drift correction | Hough Transform failure      | Adjust circle detection parameters  |

# Acknowledgments

This project is based on **TIS Imaging Control Samples** and deep-learning-based autofocus methods.

# License

This project is licensed under the MIT License. See the LICENSE file for details.

# Contact

For any issues, reach out to [**f.gielen@exeter.ac.uk**](mailto\:f.gielen@exeter.ac.uk) or open a GitHub issue.

