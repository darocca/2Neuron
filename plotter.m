pkg load signal

args = argv();
filename = args{1};

data = csvread(filename)';


t = data(1, :);
I = data(2, :);
Vm = data(3, :);
N = data(4, :);
M = data(5, :);
H = data(6, :);
I_Na = data(7, :);
I_K = data(8, :);
I_Leak = data(9, :);
postI = data(10, :);
postVm = data(11, :);

figure(2);
plot(t, Vm);
ylim ([-100, 100]);
legend ("Vm");

figure(3);
plot(t, N, 'r');
hold on
plot(t, M, 'm');
hold on
plot(t, H, 'g');
ylim([0,1]);
legend ("N", "M", "H");

figure(4);
plot(t, I_Na, 'b');
hold on
plot(t, I_K, 'k');
hold on
plot(t, I_Leak, 'r');
legend ("I Na", "I K", "I Leak");

figure(5);
V4ier1 = fft(Vm);
V4ier2 = fft(postVm);
subplot (2, 1, 1)
stem(abs(real(V4ier1)), 'r', "markersize", 0);
hold on
stem(abs(imag(V4ier1)), 'b',"markersize", 0);
axis([0, 100, 0,100000]);
xlabel ("Frequency");
ylabel ("Fourier Amplitude");
title ("Presynaptic Fourier Amplitude");
hold on;
subplot (2, 1, 2)
stem(abs(real(V4ier1)), 'r',"markersize", 0);
hold on
stem(abs(imag(V4ier2)), 'b.', "markersize", 0);
axis([0, 100, 0, 100000]);
xlabel ("Frequency");
ylabel ("Fourier Amplitude");
title ("Postsynaptic Fourier Amplitude");

%find firing rate?
function fnew = trim(f)
    for i=1:columns(f);
      if (f(1,i) < 0)
        fnew(1,i) = 0;
      else
        fnew(1,i) = f(1,i);
      end
    end
end 

[peaks1, loc1] = findpeaks(trim(Vm), "MinPeakHeight", 40, "DoubleSided");
[peaks2, loc2] = findpeaks(trim(postVm), "MinPeakHeight", 40, "DoubleSided");
totalDuration = t(columns(t));
numSpikes1 = columns(loc1);
firingRate1 = numSpikes1/totalDuration; %khZ
numSpikes2 = columns(loc2);
firingRate2 = numSpikes2/totalDuration; %khZ


figure(1);
subplot(2,2,1);
plot(t, I);
ylim ([-5,30]);
legend ("I");

subplot(2,2,2);
plot(t, Vm);
ylim ([-100, 100]);
legend ("Vm");

subplot(2,2,3);
plot(t, postI, 'r');
ylim ([-5,30]);
legend ("I");

subplot(2,2,4);
plot(t, postVm, 'r');
ylim ([-100, 100]);
legend ("Vm");

printf("Firing Rate of Neuron 1: %d%% kHz\n", firingRate1);
printf("Firing Rate of Neuron 2: %d%% kHz\n", firingRate2);

msgbox("THE PLOT THICKENS!!!");

uiwait();

