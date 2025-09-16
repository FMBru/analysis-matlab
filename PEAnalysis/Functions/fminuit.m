%FMINUIT   Multidimensional nonlinear minimization by means of the MINUIT 
%	  engine. Usage is similar (but not identical) to that of Matlab Fmins.
%	  For a description of the MINUIT commands, read the MINUIT writeup
%	  or type HELP during the interactive execution of Fminuit.
% 
%  [BestPars, Errs, Chi2, ErrMatrix] = fminuit(FunName, InitGuess, ...);
%   minimizes a scalar user function named FunName, using InitGuess as initial 
%   values of the variational parameters. Here FunName is a string and 
%   InitGuess a vector. 
%   Dots ... stand for a variable number of optional extra-arguments.
%   The user function depends on one or two variables: a vector of variational 
%   parameters (same dimension as InitGuess) and, optionally, on arbitrary 
%   structure of constant data (any data type except a string), respectively. 
%   Functions of more than 2 variables are not supported. 
%   - Note -  
%   If the user function would need to be passed many matrices of constant 
%   data, such matrices should rather be arranged into a cell array or a struct 
%   (Matlab 5 or later) or a list (Scilab 2.5 or later), and passed as a single 
%   argument to the user function.
%   
%
%  [BestPars, Errs, Chi2, ErrMatrix] = fminuit(FunName, AuxFunName, InitGuess, 
%  ...);
%   (here AuxFunName is a string) behaves as in the previous case. In addition,
%   an auxiliary function named AuxFunName is called whenever a 'CALL 5' 
%   Minuit command is issued. Such a function is passed the following 
%   argument list:  
%   (BestPars, FunName, ErrMatrix),	 or 
%   (BestPars, Data, FunName, ErrMatrix)	 if the argument Data (constant data, 
%   see below) is defined. 
%   Values of BestPars and ErrMatrix are those currently returned by Minuit. 
%   Return arguments of the auxiliary function are ignored. 
%   The interface to AuxFunName is mainly intended for logging, or to draw a 
%   run-time plot during the minimization process.
%
%
%  OPTIONAL ARGUMENTS (following InitGuess, in any order)
%
%  Data	      constant data: any data type supported by Matlab (Scilab), 
%   except a string. If this argument is supplied to Fminuit, it is passed to 
%   the user function FunName as its second argument. In the case of chi-square 
%   fitting, Data usually consists of a 3-column (or 3-row) matrix of 
%   experimental data: the independent variable, measured values of the 
%   dependent variable, and error bars.
%   - Note to the case of a scalar Data - 
%   The first occurrence of a scalar argument is interpreted as ConstError 
%   (see below). Therefore, in the case that a scalar Data is required by the 
%   user function to be minimized, ConstError has to be specified as well, and 
%   Data must follow in the argument list.
%
%  ConstError   a scalar parameter. The value returned by the user function to 
%   Fminuit is internally divided by ConstError^2. In the case of chi-square
%   fitting, ConstError has the meaning of a uniform error bar assigned to
%   the experimental points. Default value is ConstError = 1.
%
%  '-a' or 'a'	    option switch.  Interactive mode.
%   If the string constant '-a' or 'a' is supplied, Fminuit prompts for Minuit 
%   commands interactively, until 'EXIT', or 'RETURN' is issued.   
%   This is the default behaviour.
% 
%  '-b' or 'b' 	    option switch.  Batch mode.
%   If the string constant '-b' or 'b' is supplied, Fminuit runs a sequence of 
%   Minuit commands in batch mode without user control. The command sequence 
%   can be defined by the (...,'-c', MinuitCommands, ...) option, or by 
%   creating a global string variable named 'commands' in the caller workspace 
%   (see below).
%   The default sequence is 'MINIMIZE; IMPROVE' .
%
%  ..., '-c', MinuitCommands, ....	option switch followed by a string 
%   argument. Flow control. 
%   If the string constant '-c' followed by a string argument MinuitCommands 
%   is passed, Fminuit executes the Minuit commands defined in MinuitCommands 
%   when in batch mode. Defining a global non-empty string variable named 
%   'commands' in the caller workspace controls the program flow in the same 
%   way. The latter feature is kept in Fminuit for backward compatibility, 
%   and might be removed in future releases. Controlling the program flow via 
%   the command line should be preferred.
%   If both the '-c' option and the global variable commands are present, the 
%   latter is overridden.
%   Minuit commands are separated by one of the characters: ';', '\', '|', 
%   '/', or '@'; for instance,  
%   	fminuit(... ,'-c', 'set par 1 3; fix 1; mini; rel 1; mini', ...)
%   - Note - 
%   The '-c' option switch implicitly implies '-b' and is always effective; 
%   the global variable commands is only effective when batch mode is 
%   explicitly set via the '-b' (or 'b') option switch.
%
%  ..., '-n', ParameterNames, ...	option switch followed by a string 
%   argument. Parameters Names.
%   If the string constant '-n' followed by a string argument ParameterNames 
%   is passed, within Minuit, the variational parameters are given the names 
%   specified in ParameterNames, instead of the default values 'par  #1', 
%   'par  #2', ... . Defining a global non-empty string variable named 
%   'parnames' in the caller workspace controls the parameter names in the 
%   same way (backward compatibility; setting name via the command line should 
%   be preferred). 
%   If both the '-n' option and the global variable parnames are present, the 
%   latter is overridden. 
%   Parameter names are separated by commas or blanks, and must be as many as 
%   the variational parameters. If the number of names does not match the
%   number of parameters, this option is neglected.
%
%  ..., '-s', StepBounds, ....	     option switch followed by a matrix.
%   Definition of parameter steps and limits.
%   A string constant '-s' in the argument list, followed by a 2-, 3-, or 
%   4-column matrix, defines new parameter steps and/or parameter limits.
%   Defining a global non-empty matrix variable named 'stepbounds' in the 
%   caller workspace controls parameter steps and/or limits in the same way
%   (backward compatibility).
%   The '-c' switch always overrides the global variable stepbounds.
%   Format of StepBounds is the following.
%     Each line is referred to a variational parameter.
%     The first column contains indices of the parameters whose steps and/or 
%     limits are to be modified.
%     - If StepBounds is a two- or four-column matrix, the second column
%       is interpreted as the new steps, i.e. StepBounds(k,2) is the step of 
%   	the StepBounds(k,1)-th parameter;
%     - If StepBounds is a three-column matrix, columns 2 and 3 are interpreted
%   	as upper and lower limits, i.e. the StepBounds(k,1)-th parameter is 
%   	constrained between StepBounds(k,2) and StepBounds(k,3).
%     - If StepBounds is a four-column matrix, columns 2 is interpreted as new 
%       steps, column 3 and 4 as parameter limits.
%   - Note - 
%   step = 0 means that the parameter is kept fixed; 
%   a negative step signifies usage of the default step definition. 	
%
%
%  RETURN VALUES
%
%  BestPars	vector (same size as InitGuess): the optimized parameters 
%   return by Fminuit. 
%
%  Errs		vector (same size as InitGuess): the linearized (parabolic)
%   errors associated to the optimized parameters BestPars. 
%
%  Chi2		scalar: value of the user function for the variational 
%   parameters BestPars, divided by ConstError^2. Usually, this is the 
%   chi-square value corresponding to the fitting parameters returned by Fminuit
%   and the experimental data.
%
%  ErrMatrix	the error matrix, proportional to the inverse Hessian matrix.
%   ErrMatrix is a (length(BestPars), length(BestPars)) matrix, corresponding
%   row-by-row, column-by-column to the best fit parameters.
%   Fixed parameters correspond to rows and column padded with zeros.
%
%
%   UNSUPPORTED
%
%  Nonlinear errors be can calculated by issuing a MINOS Minuit commands.
%  Their values however are not interfaced to Matlab (Scilab).
%
%                               
%	Copyright (C) G. Allodi	     1996 - 2002
