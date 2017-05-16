# Kuramoto solver

Project intends to fit Kuramoto to provided data. This can be either 1-to-1 oscillator or general fitting a number of oscillators to system.   

## Fit all-to-all
In this scenario the input is an array of oscillators, i.e. S(t) = (s_1(t), s_2(t), ..., s_n(t)), and model is fit to explain what are the parameters that generated such signals.

In order to run this set configure first `model.config` file. Then, run command
```$ python kuramoto-all.py```
which will compile `stan` model, save it as `model.pkl` and results of fitted model are stored in `mode.config`.
Displaying and comparing results can be done through
```$ python plot.py```

## Fit all-to-one
The input is a linear superposition of many oscillators, thus S(T) = s_1(t)+s_2(t)+...+s_n(t).

In order to run this set configure first `model.config` file. Then, run command
```$ python kuramoto-all.py```
which will compile `stan` model, save it as `model.pkl` and results of fitted model are stored in `mode.config`.
Displaying and comparing results can be done through
```$ python plot.py```
