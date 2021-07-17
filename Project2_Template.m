% Name: Francesca Seopa
% MEC4127F: Project 2, Date 18th June 2021


pose = [0 0 0]';    %initial pose in frame {W}
dT = .01;           %sample time in seconds
T_end = 5;          %simulation duration
T = linspace(0,T_end,T_end/dT); %linearly spaced array of time points, from 0 to 5 seconds, in steps of 0.01.

p_des = [10 10]';   %desired x-y position in frame {W}. 
psi_ref = -45*pi/180; %desired final orientation for Q4. You can change this if you are testing other orientations.

cnt=0;

figure(1),clf
for ind=1:length(T)%loop makes a timestep of dT, from 0 to 5 seconds
    
    t = T(ind); %current time
    
    %ADD CONTROLLER AND POSE UPDATE FUNCTIONS HERE
    
    cnt = cnt + 1;
    if(cnt==10) %the simulation runs at 100Hz (dT=0.01), but we will plot at 10Hz
        cnt = 0;     
        
        plotTrVec = [pose(1:2); 0]; %3x1 translation vector that is made up of the current position: pose(1:2)
        plotRot = [cos(pose(3)/2) 0 0 sin(pose(3)/2)];  %construct quaternion from orientation: pose(3)
    
        xlim([0 20]),ylim([-2 20]) %fixes the x-y axes ranges
        %plots the body-frame coordinate system and also adds a visualisation of our mobile robot
        plotTransforms(plotTrVec', plotRot, "MeshFilePath", "groundvehicle.stl", "Parent", gca, "View","2D", "FrameSize", 1);
        
        %Visual indication that the controller meets the specs in Q2 (car will turn red). Comment this out when you move to Q4.
        if( norm(pose(1:2)-p_des)>1e-2 ) 
            light;  %car icon is set to white until tolerances are met. 
        end
         %Visual indication that the controller meets the specs in Q4 (car will turn red). Uncomment this when you move to Q4.
%         if( (norm(pose(1:2)-p_des)>1e-2) || (norm(psi_ref-pose(3))>5*pi/180 ) )
%             light;    %car icon is set to white until tolerances are met.
%         end
        hold on
        plot(p_des(1,end),p_des(2,end),'or','markerSize',25)    %plots our target destination as a circle.
        
        %Uncomment the line below when performing Q4. It will display the desired pose of frame {B} wrt frame {W}.
%         plotTransforms([p_des(:,end)' 0], [cos(psi_ref/2) 0 0 sin(psi_ref/2)], "Parent", gca, "View","2D", "FrameSize", 1)
        drawnow     %forces MATLAB to render the figures immediately.
    end   
    
end

// pose = [0 0 0]';    %initial pose in frame {W}
// dT = .01;           %sample time in seconds
// T_end = 5;          %simulation duration
// T = linspace(0,T_end,T_end/dT); %linearly spaced array of time points, from 0 to 5 seconds, in steps of 0.01.

// p_des = [10 10]';   %desired x-y position in frame {W}. 
// psi_ref = -45*pi/180; %desired final orientation for Q4. You can change this if you are testing other orientations.

// cnt=0;
// vl = 5;
// vr = 5;
// phi = pi/4;

// figure(1),clf
// for ind=1:length(T)%loop makes a timestep of dT, from 0 to 5 seconds
    
//     t = T(ind); %current time
    
//     %ADD CONTROLLER AND POSE UPDATE FUNCTIONS HERE
//     [vlinx, vliny, angular_phi] = forward_kinematic(vl, vr, phi);
//     %pose(1) = vlinx;
//     %pose(2) = vliny;
//     %pose(3) = angular_phi;
//     lin_vel = [vlinx; vliny; angular_phi];
//     disp("line vector");
//     disp(lin_vel);
//     rot = [cos(angular_phi) -sin(angular_phi) 0; sin(angular_phi) cos(angular_phi) 0; 0 0 1]*lin_vel;
//     disp(rot);
//     disp(pose);
//     cnt = cnt + 1;
 
//     if(cnt==10) %the simulation runs at 100Hz (dT=0.01), but we will plot at 10Hz
//         cnt = 0;     
//         %pose = pose + rot*dT;
//         plotTrVec = [pose(1:2); 0]; %3x1 translation vector that is made up of the current position: pose(1:2)
//         plotRot = [cos(pose(3)/2) 0 0 sin(pose(3)/2)];  %construct quaternion from orientation: pose(3)
    
//         xlim([0 20]),ylim([-2 20]) %fixes the x-y axes ranges
//         %plots the body-frame coordinate system and also adds a visualisation of our mobile robot
//         plotTransforms(plotTrVec', plotRot, "MeshFilePath", "groundvehicle.stl", "Parent", gca, "View","2D", "FrameSize", 1);
        
//         %Visual indication that the controller meets the specs in Q2 (car will turn red). Comment this out when you move to Q4.
//         if( norm(pose(1:2)-p_des)>1e-2 ) 
//             light;  %car icon is set to white until tolerances are met. 
//         end
//          %Visual indication that the controller meets the specs in Q4 (car will turn red). Uncomment this when you move to Q4.
// %         if( (norm(pose(1:2)-p_des)>1e-2) || (norm(psi_ref-pose(3))>5*pi/180 ) )
// %             light;    %car icon is set to white until tolerances are met.
// %         end
//         hold on
//         plot(p_des(1,end),p_des(2,end),'or','markerSize',25)    %plots our target destination as a circle.
        
//         %Uncomment the line below when performing Q4. It will display the desired pose of frame {B} wrt frame {W}.
// %         plotTransforms([p_des(:,end)' 0], [cos(psi_ref/2) 0 0 sin(psi_ref/2)], "Parent", gca, "View","2D", "FrameSize", 1)
//         drawnow     %forces MATLAB to render the figures immediately.
        
//     end   

// end

// function [vlinx, vliny, angular_phi] = forward_kinematic(vleft, vright, phi)
//     wheels_distance = 50;
//     wheel_radius = 30;
//     vlinx = (vleft + vright)*cos(phi)*wheel_radius*0.5;
//     vliny = (vleft + vright)*sin(phi)*wheel_radius*0.5;
//     angular_phi = (vleft - vright)*(wheel_radius/wheels_distance);
//     %linear_velocity = [vlinx; vliny; angular_phi];
//     %updated_pose = linear_velocity;
//     %phi = linear_velocity(3);
//     %rotation_matrix = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*linear_velocity;
//     %updated_pose = rotation_matrix*time;
// end

// % function [wheel_velocity, error_rot_vel] = position_controller(position)
// %     k1 = 10;
// %     k2 = 10;
// %     final_position = [10 10 10];
// %     error_position = final_position - position;
// %     wheel_velocity = k1*(cos(error_position(3))*error_position(1)+sin(error_position(3))*error_position(2));
// %     error_rot_vel =  k2*(arctan(error_position(1)/error_position(2)) - error_position(3));
// % end