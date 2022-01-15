clc;
close all;
clear;

R_c = 500;             %radius of the cellular cell


%BASE STATION COORDINATES IN ij COORDINATE SYSTEM
vertices_x = [R_c R_c/2 -R_c/2 -R_c -R_c/2 R_c/2 R_c];    %coordinates in x-axis
vertices_y = [0 cosd(30)*R_c cosd(30)*R_c 0 -cosd(30)*R_c -cosd(30)*R_c 0];           %coordinates in y-axis


%Coordinates of BS of first tier
vertices_BS = [1.5*R_c cosd(30)*R_c; 0 2*cosd(30)*R_c; -1.5*R_c cosd(30)*R_c];
vertices_BS = [vertices_BS; -vertices_BS];

gamma = 3.5;        %path-loss exponent
j = 1;
%%
for M = 10 : 10 : 100             %variation of number of mobile stations
    interference_users = 6*M;      %total number of users causing the interference
    
    P = -60; %Power received from any MS at its serving BS
    A = sqrt(10^((P - 30)/10));
    noise = sqrt(10^((-90 - 30)/10));       %noise power

    total_bits = 1000;       %no of iterations
    length_randm_seq = 512;      %length of sequence generated

    tx_bit = randi([0 1],M,total_bits);    %bits to be transmitted by each MS
    PN_seq = randi([0 1],M,length_randm_seq);  %PN sequence to be used by each MS

    PN = 2*PN_seq-1;
    PN_seq_bits = repmat(PN,1,total_bits);
    clear PN;
    tx_data_bits = 2*tx_bit - 1;
    tx_data_bits = repelem(tx_data_bits,1,length_randm_seq);
    tx_signal = A*tx_data_bits.*PN_seq_bits;
    clear tx_data_bits;

    net_noise = noise*randn(M,total_bits*length_randm_seq);       %overall received noise

    final_signal_sent = tx_signal + net_noise;   %final sent signal by MS
    clear net_noise tx_signal;
    rx_signal = sum(final_signal_sent);         %final received signal by BS
    users = 0;
    interference_pow = 0;    %interefering power by users from other BSs
    while(users < interference_users)
        x = -R_c + 2*R_c*rand;
        y = -R_c + 2*R_c*rand; 
        while(~inpolygon(x,y,vertices_x,vertices_y))
            x = -R_c + 2*R_c*rand;
            y = -R_c + 2*R_c*rand;
        end
        d = sqrt(x^2+y^2);        %user distance from its own BS
        users = users + 1
        BS_indices = ceil(users/M);           %The BS belongs to the user
        x = x + vertices_BS(BS_indices,1);
        y = y + vertices_BS(BS_indices,2);
        D = sqrt(x^2+y^2);           %distance of user from the reference cell
        user_PN_seq = 2*(randi([0 1],1,length_randm_seq))-1;    %PN sequence for the user
        user_bit_seq = 2*(randi([0 1],1,total_bits))-1;     %data sent by the user to the BS
        user_PN = repelem(user_PN_seq,1,total_bits);
        user_bits = repelem(user_bit_seq,1,length_randm_seq);
        user_power = A^2*d^gamma;           %Power with which the user sends the signal for perfect power control at its own BS
        rx_user_pow = user_power/D^gamma;        %power of the user received at the reference cell
        interference_pow = interference_pow + rx_user_pow;
        amp_user = sqrt(rx_user_pow);
        noise_user = noise*randn(1,total_bits*length_randm_seq);
        signal_user = amp_user*user_bits.*user_PN + noise_user;
        rx_signal = rx_signal + signal_user;
    end
    
    clear final_signal_sent;
    for m = 1 : M
        decoded_signal(m,:) = rx_signal.*PN_seq_bits(m,:);    %Multiplying with the users PN sequence to extract the data from particular user
        for i = 1 : total_bits
            rx_bits(m,i) = sign(sum(decoded_signal(m,(i-1)*length_randm_seq+1:i*length_randm_seq)));
        end
    end
    clear rx_signal decoded_signal PN_seq_bits;
    rx_bits = (rx_bits+1)/2;      %number of received bits
    error_bits = xor(rx_bits,tx_bit);    %number of error bits
    prob_of_error(j) = sum(error_bits(:))/(M*total_bits)      %error probability
    clear err_bits;
    SINR(j)=pow2db(A^2/((M-1)*A^2+interference_pow+7*M*noise^2));
    j = j + 1;
end
%%
M = 10 : 10 : 100;
figure;plot(M,prob_of_error,'Linewidth',2);title('BER vs M (Intereference from 1st Tier)');xlabel('M (no. of users)');ylabel('BER');plottools;
figure;semilogy(M,prob_of_error,'Linewidth',2);title('BER vs M (Intereference from 1st Tier)');xlabel('M (no. of users)');ylabel('BER');plottools;
figure;plot(M,SINR,'Linewidth',2);title('SINR vs M (Intereference from 1st Tier)');xlabel('M (no. of users)');ylabel('SINR in dB');plottools;