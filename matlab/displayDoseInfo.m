function [] = displayDoseInfo(filename)

load(filename);
calcBloodFracs(blood_clinical);
calcPercentKill(blood_clinical);

end
