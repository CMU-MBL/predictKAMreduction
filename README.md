# Predicting Knee Adduction Moment Response to Gait Retraining with Minimal Clinical Data
[Link to Preprint Here]

### About:
Although foot progression angle gait retraining is overall beneficial as a conservative intervention for knee osteoarthritis, knee adduction moment (KAM) reductions are not consistent across patients. Moreover, customized gait interventions are time-consuming and require instrumentation not commonly available in the clinic. We present a model that uses minimal clinical data to predict the extent of first peak KAM reduction after toe-in gait retraining. Given the lack of large public datasets that contain different gaits for the same patient, we present a method to generate toe-in gait data synthetically, and share the resultant trained model.

![Description of Data and Methods](Figures/Fig2_processDiagram.png?raw=true)

## Getting started:
Data are freely available for download at: https://simtk.org/projects/predict-kam  
With the exception of a functional data analysis component for smoothing toe-in patterns (Python), all data and code runs in MATLAB.

Scripts for processing each institution's data are contained in the folders by their name.  
The trained predictive model and learned toe-in patterns to generate synthetic toe-in gait are shared in the Models folder. 
