function [] = getTotalMU(RT_Plan_info)

mu = 0;
beams = 0;
beamNames = fieldnames(RT_Plan_info.FractionGroupSequence.Item_1.ReferencedBeamSequence);
for i = 1:numel(beamNames)
    frac = getfield(RT_Plan_info.FractionGroupSequence.Item_1.ReferencedBeamSequence, beamNames{i});
    if frac.BeamMeterset > 0
        mu = mu + frac.BeamMeterset;
        beams = beams + 1;
        
    end
end

fprintf('Total MU: %d\n', mu);
fprintf('Active Beams:%d\n', beams);

beam_on = mu / beams / (2400) * 60;
fprintf('Time on Per Beam (arc-2400):%f\n', beam_on);

beam_on = mu / beams / (1400) * 60;
fprintf('Time on Per Beam (arc-1400):%f\n', beam_on);

beam_on = mu / beams / (600) * 60;
fprintf('Time on Per Beam (static-600):%f\n', beam_on);

end
