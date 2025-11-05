# FireCA Model Repository
This repository contains multiple versions of a cellular automaton fire spread model and supporting experimental frameworks developed by different contributors. The model is informed by fire growth predictions created with Topofire which is a topographically resolved wildfire danger and drought monitoring system developed by Zachary Holden. Each folder represents a distinct phase in the model’s evolution — from early implementation to later refinements and controlled experimental testing.

---
## Max_Model/
This folder contains the code as it was when Max left the project. Max inherited the a version of the model from Adam and expanded and/or modified it. Some of the underlying logic and calculations likely originated in Adam’s implementation, though the exact portions are unclear. Comments were added in places, but certain calculations and assumptions lack complete documentation or clear reasoning.

---
## Lindsay_Model/
This folder contains the code as it was when Lindsay left the project. It builds upon and modifies the version inherited from Max. It includes detailed comments that clarify reasoning and workflow. 

---
## Toy_Model/
This folder contains a simplified experimental version of the Lindsay's CellAuto fire spread model. It was developed to test model behavior under controlled conditions using artificial (“toy”) data — for example, wind-only or slope-only scenarios — before applying the model to real fires.

---
## Computing Environment
These models are designed to be run on [Alpheus](https://alpheus_vnc7.mbdev.net/), a remote VNC-accessible virtual machine  used for model development and testing.

The code is meant to be run from the directory:
```/mnt/DataDrive1/data/zholden/VIIRS/CA_model_test/code/fireCA ```
