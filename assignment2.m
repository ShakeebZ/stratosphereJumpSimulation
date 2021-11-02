%% ENSC180-Assignment2

% Student Name 1: Shakeeb Arsalan Zacky

% Student 1 #: 301422926

% Student 1 userid (email): saz5@sfu.ca

% Student Name 2: student2

% Student 2 #: 123456782

% Student 2 userid (email): stu2 (stu2@sfu.ca)

% Below, edit to list any people who helped you with the assignment, 
%      or put ‘none’ if nobody helped (the two of) you.

% Helpers: _everybody helped us/me with the assignment (list names or put ‘none’)__

%% Instructions:
% * Put your name(s), student number(s), userid(s) in the above section.
% * Edit the "Helpers" line.  
% * Your group name should be "A2_<userid1>_<userid2>" (eg. A2_stu1_stu2)
% * Form a group 
%   as described at:  https://courses.cs.sfu.ca/docs/students
% * Replace "% <place your work here>" below, or similar, with your own answers and work.
% * Nagvigate to the "PUBLISH" tab (located on top of the editor)
%   * Click on the "Publish" dropdown and choose pdf as 
%     "Output file format" under "Edit Publishing Options..."
%   * Click "Publish" button. Ensure a report is automatically generated
% * You will submit THIS file (assignment2.m),    
%   and the PDF report (assignment2.pdf).
% Craig Scratchley, Spring 2021

%% main

function main

clf

% constants -- you can put constants for the program here
km_m_conversion = 1/1000; %1km/1000m
s_hr_conversion = 1/3600; %1hr/3600s
xlimits = [0,60];
xlimits2 = [0, 270];
xlimits3 = [255,270];
xlimits4 = [0 524.49];
g_Constant = 6.67408*10^-11;
earthMass = 5.9722*10^24; %kg
earthRadius = 6378000; %m
massOfEquipment = 27;%kg (equipment total weight)
massOf1Parachute = massOfEquipment/3; %(weight of 1 parachute)
massOfOrionParachute = 140.614;%kg
diameterOfOrionParachute = 35.3568; %m
radiusOfOrionParachute = diameterOfOrionParachute/2;%m
DeployAltitude = 2914;

Acceleration_Measured = zeros();

% prepare the data
% <place your work here>
% ... = xlsread()/csvread()/readtable()
RetrievedData = xlsread('data_clean_more_fixed_simpler');
Time_Values = RetrievedData(:,1);%seconds
Altitude = RetrievedData(:,2); %metres

Velocity_raw = RetrievedData(:,3); %kilometer/hour
Velocity_converted = (-1)*(Velocity_raw)*(s_hr_conversion)*(km_m_conversion).^-1;

Acceleration_Measured(1) = ((Velocity_converted(2))^2-(Velocity_converted(1))^2)/(2*(Altitude(2)-Altitude(1)));
for i=2:numel(Velocity_converted)
    Acceleration_Measured(i,:) = ((Velocity_converted(i))-(Velocity_converted(i-1)))/( (Time_Values(i)) - (Time_Values(i-1)));
end
Acceleration_Measured = Acceleration_Measured(:,1); %converts to positive values
Acceleration_Smooth = smoothdata(Acceleration_Measured);

%RetrievedDataUpdated = [Altitude(:,1), Velocity_converted(:,1),Time_Values(:,1)]';



% myVector(isnan(myVector))=[];
% <put here any conversions that are necessary>

%% Part 1
% Answer some questions here in these comments...
% How accurate is the model for the first portion of the minute? 

    %Our inaccuracy at the beginning of the plot was 
    %((39689.4)-(38959))/(39689.4)*100 = 1.84%
    %While the inaccuracy at 30 seconds was
    % (34407.7-34450)/34407.7*100 = -0.12%
    %The second major inaccuracy occured at approximately 29 seconds where the 
    %inaccuracy was
    %(34955 - 34845.3)/34955*100 = -0.31%
    %These inaccuracies are calculated by 
    %(MeasuredValue - ModelledValue)/MeasuredValue*100 = Inaccuracy
    %The first inaccuracy of 1.84% is due to the different starting value given
    %by the incomplete data compared to what the modelled data assumes the
    %start point to be.


% How accurate is the model for the last portion of that first minute? 
    
    %The last portion of the first minute experiences a large inaccuracy at
    %the end given by the following equation...
    %(24102-21316.8)/24102*100 = 11.56%
    %This error comes as a progressive increase from the -0.31% as the
    %model starts to diverge from the measurements over time.

% Comment on the acceleration calculated from the measured data. 
% Is there any way to smooth the acceleration calculated from the data?

    %MATLAB's Signal Processing Toolbox contains tools to make data easier
    %to interprete and read. One of which is smoothdata(), which we used so
    %that our acceleration would look a lot nicer and be less rigid. It
    %causes data to look a lot more progressive and natural by making a few
    %assumptions about how the data will behave in between plotted points.

    part = 1;
    % model Felix Baumgartner’s altitude, velocity, and acceleration for the 
    %     first minute after he jumped from 38,969.4 meters above sea level
    [T,M] = ode45(@fall, xlimits,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    VelocityModelled = AccelerationModelled(:,1);
    AccelerationModelledB = AccelerationModelled(:,2);
    
    figure(1);
    subplot(3,1,1);
    title('Part 1 - Freefall')
    plotComparisons(T, M(:,1), Altitude, Time_Values, 'Altitude (m)', 'Measured Altitude', 'Modelled Altitude')
    xlim(xlimits);
    subplot(3,1,2);
    plotComparisons(T, VelocityModelled, Velocity_converted, Time_Values, 'Velocity (m/s)', 'Measured Velocity', 'Modelled Velocity')
    xlim(xlimits);
    subplot(3,1,3);
    plotComparisons(T, AccelerationModelledB, Acceleration_Smooth, Time_Values, 'Acceleration (m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits);

%% Part 2
% Answer some questions here in these comments...
% Estimate your uncertainty in the mass that you have chosen (at the 
%    beginning of the jump). 

    % The uncertainty of the mass chosen: 114kg would be +/- 4kg, the reason
    % for this uncertainty is because we found two masses from two different
    % websites, we then took the average of the two, thus the +/- 4kg.


% How sensitive is the velocity and altitude reached after 60 seconds to 
%    changes in the chosen mass?

    % Mass changes the acceleration due to air resistance, as F/m = a, due
    %to acceleration being directly affected by mass, we can say that
    %velocity is directly affected by mass as well since we know that
    %acceleration is the time rate of change of mass.

     part = 2;      
    
     [T,M] = ode45(@fall, xlimits,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    VelocityModelled = AccelerationModelled(:,1);
    AccelerationModelledB = AccelerationModelled(:,2);
    
    figure(2);
        subplot(3,1,1);
        title('Part 2 - Air Resistance');
        plotComparisons(T, M(:,1), Altitude, Time_Values, 'Altitude (m)', 'Measured Altitude', 'Modelled Altitude')
    xlim(xlimits);
    subplot(3,1,2);
    plotComparisons(T, VelocityModelled, Velocity_converted, Time_Values, 'Velocity (m/s)', 'Measured Velocity', 'Modelled Velocity')
    xlim(xlimits);
    subplot(3,1,3);
    plotComparisons(T, AccelerationModelledB, Acceleration_Smooth, Time_Values, 'Acceleration (m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits);

%% Part 3
% Answer some questions here in these comments...
% Felix was wearing a pressure suit and carrying oxygen. Why? 
%     What can we say about the density of air in the stratosphere?
%     How is the density of air different at around 39,000 meters than it 
%     is on the ground?
	
    % If Felix were to not wear a suit, he would not survive the fall as
    % he experiences a significant g force as he opens his parachute, he
    % also needs the suit to survive the temperature conditions of the
    % fall. The oxygen tank is necessary as well due to the fact that the
    % air at that altitude is present in such small quantities that there
    % would not be enough to breathe. The density of air is much larger on
    % the ground than at 39,000 meters.

% What are the factors involved in calculating the density of air? 
%     How do those factors change when we end up at the ground but start 
%     at the stratosphere?  Please explain how calculating air density up 
%     to the stratosphere is more complicated than say just in the troposphere.

    % Pressure, Temperature, and Altitude all go into the calculation of
    % the density of rho, which the stdatmo function calculates for us
    % (Thank you!). As you go higher up into the atmosphere, there is less
    % gas available in order for us to be able to calculate the air
    % density, we must also not forget to factor pressure and temperature
    % into our calculations!

% What method(s) can we employ to estimate [the ACd] product? 

    % We thought of using the formula: CDA = 2*m*(g - a*z)/(rho*v^2), however,
    % whenever we implemented this into our code, it fails to work when the
    % ODE45 function calls the fall function and accesses the drag. We
    % discovered that the velocity retrived would (for some reason) equal
    % to zero, when using the ODE45 call. This resulted in us get a CdA of
    % infinity or undefined, and so we decided to find constants for the
    % CdA instead. We first assumed him to be a sphere and so the Cd of 0.5
    % was found. We then needed to find his area, we searched around and
    % found another link with a calculated area in it (0.8m^2). Using these
    % values we calculated CdA by multiplying Cd and A, giving us 0.4CdA.
    %link to Cd constants: https://www.engineeringtoolbox.com/drag-coefficient-d_627.html
    %link to Area of body: https://physicstoday.scitation.org/doi/full/10.1063/PT.3.2357#:~:text=Here%20m%20is%20Baumgartner%27s%20118,the%20altitude%2Ddependent%20gravitational%20acceleration

% What is your estimated [ACd] product?
    
    % Our estimated CdA product is: 0.4ACd.

% [Given what we are told in the textbook about the simple drag constant, b,] 
%   does the estimate for ACd seem reasonable?
    
    % Considering it is less than 1, and near the given drag constant of 0.2,
    % we would assume this estimate seems reasonable.

part = 3;
     [T,M] = ode45(@fall, xlimits2,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    VelocityModelled = AccelerationModelled(:,1);
    AccelerationModelledB = AccelerationModelled(:,2);
    
    figure(3);
        subplot(3,1,1);
        title('Part 3 - Drag');
        plotComparisons(T, M(:,1), Altitude, Time_Values, 'Altitude (m)', 'Measured Altitude', 'Modelled Altitude')
    xlim(xlimits2);
    subplot(3,1,2);
    plotComparisons(T, VelocityModelled, Velocity_converted, Time_Values, 'Velocity (m/s)', 'Measured Velocity', 'Modelled Velocity')
    xlim(xlimits2);
    subplot(3,1,3);
    plotComparisons(T, AccelerationModelledB, Acceleration_Smooth, Time_Values, 'Acceleration (m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits2);
%% Part 4
% Answer some questions here in these comments...
% What is the actual gravitational field strength around 39,000 meters? 
    
    %The actual gravitational field strength around 39,000 meters is
    %-9.673777

% How sensitive is the altitude reached after 4.5 minutes to simpler and 
%   more complicated ways of modelling the gravitational field strength? 
% if we were to use -9.81m/s^2 as our acceleration and we were to compare it to
% the values if we had used the value from the gravitational field strength
% formula (-9.6737777m/s^2), the difference can be calculated using:
% 1/2at^2 if we were to only worry about the acceleration. Using this, the
% displacement of using -9.81m/s^2 as our acceleration would give:
% 357,574.5m after 4.5minutes, whereas if we were to us the changing
% gravitational field we would have gotten a displacment of: 325.606.2m.
% This is a difference of 4,965.3 meters. 


% What other changes could we make to our model? Refer to, or at least 
%   attempt to explain, the physics behind any changes that you propose. 

    %The change in area as Felix opens his parachute and correspondingly
    %the changing value of CdA. As area increases the drag increases which
    %changes our acceleration, and thus our velocity graphs. This would not
    %drastically change our graphs but there would be a slight change over
    %time.

% What is a change that we could make to our model that would result in 
%   insignificant changes to the altitude reached after 4.5 minutes? 
    
    % Even if we change the mass of Felix or his equipment, this would
    % result in no change to the model. This is due to the fact that
    % calculating the gravitational strength is not affected by the mass of
    % the secondary object/person compared to the larger mass. In our case:
    % the Earth.

% How can we decide what change is significant and what change is 
%   insignificant?

    % We can determine what is significant/insignificant by looking at how
    % the values affect our calculations. For example, had his velocity
    % increase or decreased, this would have drastically altered the drag
    % force. As well, if we has chosen him to be represented as another
    % object instead of the sphere, his Cd would differ.

% [What changes did you try out to improve the model?  (Show us your changes
%   even if they didn't make the improvement you hoped for.)]
    
    %We modelled CdA as two different values, one of which was due to Felix
    %and the other due to Felix and the parachute being open. This change
    %can be seen in the graphs and in our drag function.



part = 4;
[T,M] = ode45(@fall, xlimits2,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    VelocityModelled = AccelerationModelled(:,1);
    AccelerationModelledB = AccelerationModelled(:,2);

figure(4);
    subplot(3,1,1);
    title('Part 4 - Gravity');
    plotComparisons(T, M(:,1), Altitude, Time_Values, 'Altitude (m)', 'Measured Altitude', 'Modelled Altitude')
    xlim(xlimits2);
    subplot(3,1,2);
    plotComparisons(T, VelocityModelled, Velocity_converted, Time_Values, 'Velocity (m/s)', 'Measured Velocity', 'Modelled Velocity')
    xlim(xlimits2);
    subplot(3,1,3);
    plotComparisons(T, AccelerationModelledB, Acceleration_Smooth, Time_Values, 'Acceleration (m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits2);

%% Part 5
% Answer some questions here in these comments...
% At what altitude does Felix pull the ripcord to deploy his parachute? 

    %Felix deploys the parachute at approximately 2918 meters +- 100 meters
    %as seen from the video. The parachute is fully deployed around 2500
    %meters.

% Recalculate the CdA product with the parachute open, and modify your 
%   code so that you use one CdA product before and one after this altitude. 
%   According to this version of the model, what is the maximum magnitude 
%   of acceleration that Felix experiences? 

    % According our model, the maximum magnitude Felix experiences is:
    % 101.253 m/s^2

%   How safe or unsafe would such an acceleration be for Felix?

    % This acceleration is much greater than what any human experiences day by
    % day, This is equivelent to 10.32G's, this is 1G greater than that of what
    % fight pilots experience. Fighter pilots experience up to 9G's, Although
    % greater than what the pilots expereinced, it is a safe acceleration of
    % Felix to feel, as around 18G's is when acceleration becomes lethal. -->
    % according to this link: https://www.pbs.org/wgbh/nova/warplanes/gforces.html#:~:text=With%20the%20development%20of%20faster,of%20gravity%20at%20sea%20level.

part = 5;
%Make a single acceleration-plot figure that includes, for each of the 
%model and the acceleration calculated from measurements, the moment when 
%the parachute opens and the following 10 or so seconds. If you have 
%trouble solving this version of the model, just plot the acceleration 
%calculated from measurements. 
    [T,M] = ode45(@fall, xlimits4,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    AccelerationPart5 = AccelerationModelled(:,2);
    figure(5)
    title('Part 5 - Parachute Deployment');
    plotComparisons(T, AccelerationPart5, Acceleration_Smooth, Time_Values, 'Acceleration (m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits3);

%% Part 6 
% Answer some questions here in these comments...
% How long does it take for Felix’s parachute to open?

    %It takes approximately 7.974 seconds for Felix's parachute to open as
    %calculated from the Redbull Stratos Jump Video.

part = 6;

%Redraw the acceleration figure from the previous Part but using the new 
%   model. Also, using your plotting function from Part 1, plot the 
%   measured/calculated data and the model for the entire jump from 
%   stratosphere to ground.
% <place your work here>

[T,M] = ode45(@fall, xlimits4,[38969.4,0]); 
    AccelerationModelled = zeros(2, numel(T));
    for i=1:numel(T)
        AccelerationModelled(1:2,i) = fall(T(i),M(i,:));
    end
    AccelerationModelled = AccelerationModelled';
    VelocityModelled = AccelerationModelled(:,1);
    AccelerationModelledB = AccelerationModelled(:,2);
    
    figure(6)
    subplot(4, 1, 1);
    title('Part 6 - Results');
    plotComparisons(T, M(:,1), Altitude, Time_Values, '(m)', 'Measured Altitude', 'Modelled Altitude')
    xlim(xlimits4);
    subplot(4,1,2);
    plotComparisons(T, VelocityModelled, Velocity_converted, Time_Values, '(m/s)', 'Measured Velocity', 'Modelled Velocity')
    xlim(xlimits4);
    subplot(4,1,3);
    plotComparisons(T, smoothdata(AccelerationModelledB), Acceleration_Smooth, Time_Values, '(m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits4);
    subplot(4, 1, 4);
    plotComparisons(T, smoothdata(AccelerationModelledB), Acceleration_Smooth, Time_Values, '(m/s^2)', 'Measured Acceleration', 'Modelled Acceeleration')
    xlim(xlimits3);
    
%% nested functions  
% nested functions below are required for the assignment.  
% see Appendix B of Physical Modeling in MATLAB for discussion of nested functions

function res = fall(t, X)
    %FALL Function takes in parameters of t and X, function returns
    %velocity and acceleration.
    %   
    %when called, fall takes in two parameters, t, and X. X is then
    %seperated into X(1) and X(2) where X(1) is the first row of the vector
    %and X(2) is the second row. The row vectors are entered into functions
    %relating to acceleration(t, y, v), gravityEst(y), mass(t, v) and drag(t, y, v, m).
    %Where the value of acceleration determined through a series of
    %calculations.

    % do not modify this function unless required by you for some reason! 

    y = X(1); % the first element is position
    v = X(2); % the second element is velocity

    dydt = v; % velocity: the derivative of position w.r.t. time
    dvdt = acceleration(t, y, v); % acceleration: the derivative of velocity w.r.t. time

    res = [dydt; dvdt]; % pack the results in a column vector
end

function res = acceleration(t, y, v)
    % ACCELERATION Function takes in parameters t, y, and v, and returns
    % acceleration
    % input...
    %   t: time
    %   y: altitude
    %   v: velocity
    % output...
    %   res: acceleration
    %
    %If the part is 1, then the acceleration that is returned will always
    %be the negative value of the gravityEst function.
    %If the part is any value other than 1 then the mass function will be
    %called along with the drag function. Mass takes in parameters t and v 
    %and returns m. Drag takes in parameters t, y, v, and m. It returns b
    %which the acceleration function uses to calculate our value.

    % do not modify this function unless required by you for some reason! 

    grav = gravityEst(y); 

    if part == 1 % variable part is from workspace of function main.
        res = -grav;
    else
        m = mass(t, v);
        b = drag(t, y, v, m);

        f_drag = -b * v^2 * sign(v);
        a_drag = f_drag / m;
        res = -grav + a_drag;
    end
end

% Please paste in or type in code into the below functions as may be needed.

function grav = gravityEst(y)
    % estimate the acceleration due to gravity as a function of altitude, y
    g_SEA = 9.807;  % gravity at sea level in m/s^2

    if part <= 3
        grav = g_SEA;
    else
        grav = (((g_Constant)*(earthMass)) / (earthRadius+y).^2);
        
     end
end

function res = mass(t, v)
    mass_of_person = 114.0; % (kg) this value was determined by summing the weight of Felix, and the weight of his equipment
    net_weight = mass_of_person; 
    
    res = net_weight;  
end

function res = drag(t, y, v, m)
    %DRAG Function takes in parameters t, y, v, and m.
    %If it is called with part equal to 2, then it will always return 0.2;
    %however, if part does not equal 2 then it calculates rho using the
    %stdatmo function. A constant CdA is chosen depending on the y variable
    %which is used with rho to determine the drag.
    MassRatioParachutes = massOf1Parachute/massOfOrionParachute;
    %Assuming the mass ratio will constitute the difference
    %in diameter of the parachutes: the Orian Parachute used
    %in Nasa vehicles, and the parachute Felix used.
    ParachuteRadius = MassRatioParachutes*radiusOfOrionParachute;
     if y <= DeployAltitude && part > 5
        if y >= 2850
            ParachuteRadius = (MassRatioParachutes*radiusOfOrionParachute)/(2);
        end
     end
     
    AreaOfParachute = pi*(ParachuteRadius)^2;
    
    if part == 2
        res = 0.2;
    else
        % air resistance drag = 1/2*rho*c_d*a = 1/2*rho*CdA
        if y>0
            rho = stdatmo(y);
        else
            rho = stdatmo(abs(y));
        end
        if y <= DeployAltitude && part > 3
            Cd = 1.25; %This value is the average coefficient for parachutes for 
                       %a bottled rocket and a spacecraft parachute. The
                       %reason for averaging it is to find the middle
                       %grounds between the the smallest parachute and the
                       %largest ones which nasa uses to land its
                       %spacecrafts.
                       %https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktvrecv.html#:~:text=Typical%20values%20of%20drag%20coefficient,produces%20a%20lower%20terminal%20velocity.
            CdA = Cd*(AreaOfParachute);    
        else
        AreaOfBody = 0.8; % <--This area was not calculated by us, value
                          % was found in this link:
                          % https://physicstoday.scitation.org/doi/full/10.1063/PT.3.2357#:~:text=Here%20m%20is%20Baumgartner%27s%20118,the%20altitude%2Ddependent%20gravitational%20acceleration
                          
        Cd = 0.5; % <-- assuming He is considered a sphere
                       % https://www.engineeringtoolbox.com/drag-coefficient-d_627.html
                       % ^ link to the coefficients of Cd
        CdA = Cd*AreaOfBody;
        end
        drag = (CdA.*rho);
        res = drag;
    end
end

%% Additional nested functions
% Nest any other functions below.  
%Do not put functions in other files when you submit, except you can use
%    the stdatmo function in file stdatmo.m which has been provided to you.
       function plotComparisons(TData, ModelledData, Measured, Time, Title, Legend1, Legend2)
        LegendArray = {Legend1, Legend2};
        hold on
        xlabel('Time (s)');
        plot(Time, Measured, 'r-')
        plot(TData,ModelledData,'b-'); 
        ylabel(Title);
        legend(LegendArray, 'Location', 'northeast');
       end   
% end of nested functions
end % closes function main.  