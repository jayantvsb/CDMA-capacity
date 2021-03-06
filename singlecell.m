clc;
close all;
clear;

j = 1;     %there is only one cell i.e. reference cell
for M = 10 : 10 : 100         %variation of number of mobile stations
    total_bits = 1000;            %No of iterations
    length_randm_seq = 512;          %length of sequence generated
    P = -60;               %Power received from any MS at its serving BS
    A = sqrt(10^((P - 30)/10));
    noise = sqrt(10^((-90 - 30)/10));             %noise power

    tx_bits = randi([0 1],M,total_bits);    %bits to be transmitted by each MS
    PN_seq = randi([0 1],M,length_randm_seq);  %PN sequence to be used by each MS

    PN = 2*PN_seq-1;
    PN_seq_bits = repmat(PN,1,total_bits);
    clear PN;
    tx_data_bits = 2*tx_bits - 1;
    tx_data_bits = repelem(tx_data_bits,1,length_randm_seq);
    tx_signal = A*tx_data_bits.*PN_seq_bits;
    clear tx_data_bits;

    net_noise = noise*randn(M,total_bits*length_randm_seq);  %overall noise received

    final_signal_sent = tx_signal + net_noise;   %final signal sent by the MS
    clear net_noise tx_signal;
    rx_signal = sum(final_signal_sent);        %final received signal by BS
    clear final_signal_sent;
    for m = 1 : M
        decoded_signal(m,:) = rx_signal.*PN_seq_bits(m,:);       %decoded signal onto the receiver side
        for i = 1 : total_bits
            rx_bits(m,i) = sign(sum(decoded_signal(m,(i-1)*length_randm_seq+1:i*length_randm_seq)));
        end
    end
    clear rx_signal decoded_signal PN_seq_bits;
    rx_bits = (rx_bits+1)/2;             %total number of receive bits
    error_bits = xor(rx_bits,tx_bits);    %number of error bits
    prob_of_error(j) = sum(error_bits(:))/(M*total_bits);           %error probability in the number of bits
    clear error_bits;
    SINR(j)=pow2db(A^2/((M-1)*A^2+M*noise^2));          %signal to noise ratio
    j = j + 1;
end
%%
M = 10:10:100;
figure;plot(M,prob_of_error,'Linewidth',2);title('BER vs M');xlabel('M (no. of users)');ylabel('BER');plottools;
figure;semilogy(M,prob_of_error,'Linewidth',2);title('BER vs M');xlabel('M (no. of users)');ylabel('BER');plottools;
figure;plot(M,SINR,'Linewidth',2);title('SINR vs M');xlabel('M (no. of users)');ylabel('SINR in dB');plottools;