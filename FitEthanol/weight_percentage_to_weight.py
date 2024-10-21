import math
import csv

def lerp(a, b, t):
    return (1 - t) * a + t * b

# From the wikipedia data page of ethanol
def ethanol_density_func(temperature):
    return -8.461834e-4 * temperature + 0.8063372
    
# From https://holzmann-cfd.com/community/blog-and-tools/cae-blog/thermophysical-properties-water
def water_density_func(temperature):
    temperature_kelvin = temperature + 273.15
    return -0.00365471 * temperature_kelvin ** 2 + 1.93017 * temperature_kelvin + 746.025
            
def handle_mixing_input(name, recommended):
    while True:
        try:
            return float(input(f"Please weigh {recommended} g of {name} and enter real weighed amount: "))
        except:
            print("Input a real number please")
            
def write_table_for_results(path, results):
    print(results)
    print(results[0])
    with open(path, newline='', mode='w') as csv_file:
        writer = csv.writer(csv_file, delimiter=' ')
        writer.writerow(['n', 'kg weighed'])
        
        for i in range(len(results[0])):
            writer.writerow([str(i), str(results[0][i])])
            
        writer.writerow(['', ''])
        writer.writerow(['Excess amount: ', str(results[1])])

def mix_helper(number, name, max_weight):
    steps = []
    while number > 0:
        amount = handle_mixing_input(name, min(max_weight, number * 1000)) / 1000
        number -= amount
        steps.append(amount)
        
    print(f"{abs(number) * 1000} g too much of {name}.")
    print(" --- ")
    
    return [steps, abs(number)]

# Everyting in SI
ethanol_temperature = 25
water_temperature = 25
# volume = 0.9e-3
volume = 0.000448
ethanol_density = ethanol_density_func(ethanol_temperature) * 1000
water_density = water_density_func(water_temperature)
weight_percentage = 0.12615

# This is the maximum weight of the scale in g
max_weight = 5000

weight_water = 1.0 / ((weight_percentage) / (volume * ethanol_density * (1 - weight_percentage)) + 1.0 / (water_density * volume))
weight_ethanol = (weight_percentage * weight_water) / (1 - weight_percentage)    

# Result checks out
print(f"Weight of water: {weight_water} kg")
print(f"Weight of ethanol: {weight_ethanol} kg")

result_water = mix_helper(weight_water, "water", max_weight)
actual_weight_water = weight_water + result_water[1]
actual_weight_ethanol = (weight_percentage * actual_weight_water) / (1 - weight_percentage)
result_ethanol = mix_helper(actual_weight_ethanol, "ethanol", max_weight)
actual_weight_ethanol += result_ethanol[1]

print(result_water)

actual_weight_percentage = actual_weight_ethanol / (actual_weight_ethanol + actual_weight_water)

# Outputting things
experiment_parameters_path = 'experiment_params.csv'
water_weighing_data_path = 'water_weighing_data.csv'
ethanol_weighing_data_path = 'ethanol_weighing_data.csv'

write_table_for_results(water_weighing_data_path, result_water)
write_table_for_results(ethanol_weighing_data_path, result_ethanol)

with open(experiment_parameters_path, newline='', mode='w', encoding="utf-8") as csv_file:
    writer = csv.writer(csv_file, delimiter=' ')
    writer.writerow(['T_e', str(ethanol_temperature)])
    writer.writerow(['T_w', str(water_temperature)])
    writer.writerow(['V', str(volume)])
    writer.writerow(['ρ_e', str(ethanol_density)])
    writer.writerow(['ρ_w', str(water_density)])
    writer.writerow(['ξ_theory', str(weight_percentage)])
    writer.writerow(['m_e_theory', str(weight_ethanol)])
    writer.writerow(['m_w_theory', str(weight_water)])
    writer.writerow(['ξ_reality', str(actual_weight_percentage)])
    writer.writerow(['m_e_reality', str(actual_weight_ethanol)])
    writer.writerow(['m_w_reality', str(actual_weight_water)])
    writer.writerow(['m_scale_max', str(max_weight / 1000)])

input("Program ended")