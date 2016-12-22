'''
This script contains the specific parameters to make the sub-plots pretty for each of the 4 plots
for all objects. 
'''

# define dictionary
galaxies_subplt_dict = {
                        'mrk960':{
                                  'NUV':{
                                         'plt_number' : 1,
                                         'xlims' : [1590.0, 2100.0],
                                         'lines' : [1661, 1666, 1909],
                                         'legends' : ['O III]', 'O III]', 'C III]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [32],
                                         'subycoord_list' : [68],
                                         'side_list' : ['center'],
                                         'lineymax_list' : [0.8e-15],
                                         'idx_list' : [0]
                                         },
                                  'Blue':{
                                         'plt_number' : 2,
                                         'xlims' : [3700.0, 5100.0],
                                         'lines' : [3727, 3869, 3967, 4069, 4076, 4102, 4340, 4363, 
                                                    4861, 4959, 5007],
                                         'legends' : ['[O II]', '[Ne III]', '[Ne III]', '[S II]', '[S II]', 
                                                      'Hg', 'Hd', '[O III]', 'Hb', '[O III]', '[O III]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [0, 0, 0, 0, 0, -4, 0, 0, 0, 0, 0],
                                         'subycoord_list' : [158, 128, 78, 68, 68, 48, 44, 38, 58, 58, 58],
                                         'side_list' : ['left', 'left', 'left', 'right', 'left', 'left', 'left', 
                                                        'left', 'right', 'right', 'right'],
                                         'lineymax_list' : [0.7e-15, 0.55e-15, 0.3e-15, 0.3e-15, 0.3e-15, 0.2e-15, 
                                                            0.17e-15, 0.12e-15, 0.2e-15, 0.2e-15, 0.2e-15],
                                         'idx_list' : [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]
                                         },
                                  'Red':{
                                         'plt_number' : 3,
                                         'xlims' : [6250.0, 6750.0],
                                         'lines' : [6300, 6312, 6563, 6678, 6716],
                                         'legends' : ['[O I]', '[S III]', 'Ha', 'He I', '[S II]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [32, -4, -70, 32, -48],
                                         'subycoord_list' : [60, 18, 18, 60, 18],
                                         'side_list' : ['center', 'left', 'left', 'right', 'left'],
                                         'lineymax_list' : [1.2e-16, 0.3e-16, 0.25e-16, 1.2e-16, 0.25e-16],
                                         'idx_list' : [0, 0, 1, 0, -1]
                                         },
                                  'NIR':{
                                         'plt_number' : 4,
                                         'xlims' : [9000.0, 9600.0],
                                         'lines' : [9069, 9531],
                                         'legends' : ['[S III]', '[S III]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [-18, -10],
                                         'subycoord_list' : [55, 32],
                                         'side_list' : ['left', 'center'],
                                         'lineymax_list' : [1.3e-16, 0.5e-16],
                                         'idx_list' : [-1, 0]
                                         }
                                  },
                        
                        'sbs0218':{
                                  'NUV':{
                                         'plt_number' : 1,
                                         'xlims' : [1590.0, 2100.0],
                                         'lines' : [1661, 1666, 1909],
                                         'legends' : ['O III]', 'O III]', 'C III]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [32],
                                         'subycoord_list' : [101],
                                         'side_list' : ['center'],
                                         'lineymax_list' : [0.5e-15],
                                         'idx_list' : [1, 0, 1]
                                         },
                                  'Blue':{
                                         'plt_number' : 2,
                                         'xlims' : [3700.0, 5100.0],
                                         'lines' : [3727, 3869, 3967, 4069, 4076, 4102, 4340, 4363, 
                                                    4861, 4959, 5007],
                                         'legends' : ['[O II]', '[Ne III]', '[Ne III]', '[S II]', '[S II]', 
                                                      'Hg', 'Hd', '[O III]', 'Hb', '[O III]', '[O III]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [0, 0, 0, 0, 0, -4, 0, 0, 0, 0, 0],
                                         'subycoord_list' : [100, 128, 78, 68, 68, 48, 44, 38, 58, 58, 58],
                                         'side_list' : ['left', 'left', 'left', 'right', 'left', 'left', 'left', 
                                                        'left', 'right', 'right', 'right'],
                                         'lineymax_list' : [0.8e-15, 0.65e-15, 0.3e-15, 0.3e-15, 0.3e-15, 0.2e-15, 
                                                            0.17e-15, 0.12e-15, 0.2e-15, 0.2e-15, 0.2e-15],
                                         'idx_list' : [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0]
                                         },
                                  'Red':{
                                         'plt_number' : 3,
                                         'xlims' : [6250.0, 6750.0],
                                         'lines' : [6300, 6312, 6563, 6678, 6716],
                                         'legends' : ['[O I]', '[S III]', 'Ha', 'He I', '[S II]'],
                                         'fmin' : -0.25e-15,
                                         'subxcoord_list' : [32, -4, -70, 32, -48],
                                         'subycoord_list' : [60, 18, 18, 60, 18],
                                         'side_list' : ['center', 'left', 'left', 'right', 'left'],
                                         'lineymax_list' : [1.2e-16, 0.3e-16, 0.25e-16, 1.2e-16, 0.25e-16],
                                         'idx_list' : [0, 0, 1, 0, -1]
                                         },
                                  'NIR':{
                                         'plt_number' : 4,
                                         'xlims' : [9000.0, 9600.0],
                                         'lines' : [9069, 9531],
                                         'legends' : ['[S III]', '[S III]'],
                                         'fmin' : -0.65e-15,
                                         'subxcoord_list' : [-18, -10],
                                         'subycoord_list' : [55, 32],
                                         'side_list' : ['left', 'center'],
                                         'lineymax_list' : [1.3e-16, 0.5e-16],
                                         'idx_list' : [-1, 0]
                                         }
                                  }


                        }

