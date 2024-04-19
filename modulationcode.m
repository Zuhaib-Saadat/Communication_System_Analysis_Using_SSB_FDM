 pause(1)
 
 disp('Start speaking.')
recordblocking(recObj2,5); 
% recordblocking(OBJ, T) records for length of time, T, in seconds;
%                             does not return until recording is finished.
disp('End of Recording.');
%play(recObj2);
mySpeech = getaudiodata(recObj2); % returns the recorded audio data as a double array
% getaudiodata(OBJ, DATATYPE) returns the recorded audio data in
%      the data type as requested in string DATATYPE.  Valid data types
%      are 'double', 'single', 'int16', 'uint8', and 'int8'.
%Write audio file
audiowrite('m2.wav',mySpeech,Fs);

 pause(1)
 
 disp('Start speaking.')
recordblocking(recObj3,5); 
% recordblocking(OBJ, T) records for length of time, T, in seconds;
%                             does not return until recording is finished.
disp('End of Recording.');
%play(recObj3);   %PLAYS AUDIO BACK TO YOU
mySpeech = getaudiodata(recObj3); % returns the recorded audio data as a double array
% getaudiodata(OBJ, DATATYPE) returns the recorded audio data in
%      the data type as requested in string DATATYPE.  Valid data types
%      are 'double', 'single', 'int16', 'uint8', and 'int8'.
%Write audio file
audiowrite('m3.wav',mySpeech,Fs);



