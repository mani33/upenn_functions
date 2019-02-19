function set_matlab_window_title(str)
% Sets a custom title for the matlab main window so that when you run
% multiple matlabs, you know which one to go to.
jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jDesktop.getMainFrame.setTitle(str);