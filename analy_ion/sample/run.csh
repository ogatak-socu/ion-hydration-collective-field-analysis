#
../analy_ion/cal_vector -op nacl_sample.1.in
../analy_ion/cal_vector -op nacl_sample.2.in
../analy_ion/cal_vector -op nacl_sample.3.in

../analy_ion/accumulate_vector -op nacl_sample.acc.in

../analy_ion/cal_charge -op nacl_sample.1.in
../analy_ion/cal_charge -op nacl_sample.2.in
../analy_ion/cal_charge -op nacl_sample.3.in

../analy_ion/accumulate_charge -op nacl_sample.acc.in


