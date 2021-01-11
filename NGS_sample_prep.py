from opentrons import protocol_api
import numpy as np
import math 

def get_values(*names):
    import json
    _user_input = json.loads("""{"num_samples":5,"MiSeq_Sequencer":"Geena","qPCR_quant":[500, 50, 50, 50, 50],"tapestation":[500, 0, 50, 50, 50], "elution_buffer":50, "volume":50, "coverage":[20, 20, 20, 20, 20], "diversity":[100, 100, 100, 100, 100]}""")
    return [_user_input[n] for n in names]

metadata = {
    'protocolName': 'Primer',
    'author': 'Sebastian <sebastian.lewis@labgeni.us',
    'apiLevel': '2.5'
}

def run(protocol: protocol_api.ProtocolContext):
    [num_samples, qPCR_quant, tapestation, elution_buffer, volume, diversity, coverage, MiSeq_Sequencer] = get_values(
        'num_samples', 'qPCR_quant', 'tapestation', 'elution_buffer', 'volume', 'diversity', 'coverage', 'MiSeq_Sequencer')
    
    input_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '6')
    reaction_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '4')
    reaction_reagents = protocol.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', '3')
    destination_2nd = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '5')
    normalised_pooling = protocol.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', '8')
    pooled_denatured= protocol.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', '7')
    tiprack = [
        protocol.load_labware('opentrons_96_tiprack_300ul', slot)
                for slot in ['1', '9']
                ]
    p10_tiprack = [
        protocol.load_labware('opentrons_96_tiprack_10ul', slot)
        for slot in ['2']
        ]
    
    p300 = protocol.load_instrument(
        'p300_single', 'right', tip_racks=tiprack)
    p10 = protocol.load_instrument('p10_single', 'left', tip_racks=p10_tiprack)

    empty_list = []

    for _ in range(1):
        if qPCR_quant == 0:
            average_values = [((x) for x, y in list(zip(qPCR_quant, tapestation)))]
            empty_list.append(average_values)
        elif tapestation == 0:
            average_values = [((y) for x, y in list(zip(qPCR_quant, tapestation)))]
            empty_list.append(average_values[:num_samples])
        else:
            average_values = [(((x+y)/2)) for x, y in list(zip(qPCR_quant, tapestation))]
            empty_list.extend(average_values[:num_samples])
   
    print(average_values)
   
    protocol.comment("Welcome to the NGS normalisation, pooling and denaturation protocol")
    
    average_values_array = np.array(average_values)
    ten_fold_dilution = (average_values_array/10)
    nm_to_ul = ((volume*4/average_values_array))
    dilution_buffer_to_add = (elution_buffer-nm_to_ul)
    dilution_buffer_to_adds = list(dilution_buffer_to_add)
    total_volume_combined = (dilution_buffer_to_adds+nm_to_ul)
    average_values_list=list(nm_to_ul)
    total_volume=np.array(total_volume_combined)

    for average_vales in average_values_list:
        if any(average_values) < 1:
            numpy.ceil(average_values)

    protocol.comment('Normalisng the {} samples'.format(num_samples))

# Uses the calculations above to find 4nm_to_ul of each library

    p10.transfer(
        average_values_list[:num_samples],
        input_plate.wells()[:num_samples],
        reaction_plate.wells()[:num_samples],
        new_tip='always')

# Uses the calculations above to add elution buffer to make each sample 50 uL
        
    p300.transfer(
        dilution_buffer_to_adds[:num_samples],
        reaction_reagents['A1'],
        reaction_plate.wells()[:num_samples],
        new_tip='always')

#Pooling

    protocol.comment('Currently pooling the {} samples'.format(num_samples))

#Pooling calculations

    diversitys=np.array(diversity)
    coverages=np.array(coverage)
    number_of_reads=np.multiply(coverages, diversitys)
    percentage_of_run=(((number_of_reads)/4000000)*100)
    percentage_of_run_per_library=np.multiply(total_volume, percentage_of_run)
    percentage_of_run_per_library_list=list(percentage_of_run_per_library)

#Pooling each library based on percentage of run

    p10.transfer(
        percentage_of_run_per_library_list[:num_samples],
        reaction_plate.wells()[:num_samples],
        normalised_pooling['A1'],
        new_tip='always')
        
#denaturation

# Ths makes the NaOH and dilutes to 1 nM

    protocol.comment('Currently Denaturing the {} samples'.format(num_samples))

    p300.transfer(
        200,
        reaction_reagents['A1'],
        reaction_reagents['C1'],
    )

    p300.transfer(
        800, 
        reaction_reagents['A1'],
        reaction_reagents['C1'],
    )

# Adds 5 uL of library and 5 uL of NaOH to denature
   
    p10.transfer(
        5,
        normalised_pooling['A1'],
        pooled_denatured['A1'],
        new_tip='always')
    
    p10.transfer(
        5,
        reaction_reagents['C1'],
        pooled_denatured['A1'],
        new_tip='always')
    
    protocol.pause('Vortex briefly and then centrifuge at 280 × g for 1 minute.')
    protocol.delay(minutes=5) 

# Adds 990uL of hybridisation buffer to denatured librarys
    
    p300.transfer(
        990, 
        reaction_reagents['C1'],
        pooled_denatured['A1'],
        new_tip='always')

# Dilutes the denatured library even further based on the sequencer used. 
    
    if MiSeq_Sequencer == "Geena":
        
        p300.transfer(
        165, 
        input_plate.wells()[:num_samples],
        destination_2nd.wells()[:num_samples],
        new_tip='always')
        
        p300.transfer(
        435, 
        input_plate.wells()[:num_samples],
        destination_2nd.wells()[:num_samples],
        new_tip='always')
    
    else:

        p300.transfer(
        225, 
        input_plate.wells()[:num_samples],
        destination_2nd.wells()[:num_samples],
        new_tip='always')
        
        p300.transfer(
        375, 
        input_plate.wells()[:num_samples],
        destination_2nd.wells()[:num_samples],
        new_tip='always')
