function [] = stepResp(step_offset,step_Amp,Tf,stepT,scaleY,input_index,output_index,A,B,C,D)

sys               = ss(A,B,C,D);
figure
opt               = stepDataOptions;
opt.StepAmplitude = step_Amp;
t                 = 0:stepT:Tf-step_offset;
y                 = step(sys,t,opt);
if step_offset == 0
    t     = 0:stepT:Tf;
    y     = step(sys,t,opt);
    plot(t, y(:, output_index, input_index));   %Remember: you have n_outputs and m_inputs such that size(D_tot0) = [n,m] where m -> uk_s (the global inputs)
    title(sprintf('output (%s) after a step change of input (%s) of magnitude %.2f', output_index, input_index, opt.StepAmplitude));
    xlabel('Time');
    ylabel(sprintf('%s', output_index));
else
    t1   = 0:stepT:(step_offset-stepT);
    t2   = step_offset:stepT:Tf;
    t    = [t1 t2];
    y1   = scaleY*ones(length(t1),1);
    y2   = [y1; scaleY+y(:, output_index, input_index).*scaleY];
    plot(t,y2)
end

end