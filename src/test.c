#include "funs.h"
#include "heap.h"
#include "gmc.h"

#include <stdio.h>
#include <sys/time.h>

//static double v1[400] = { /*First row*/ -1.2887, -1.3371, -1.3821, -1.2948, -1.5536, -1.5409, -1.5559, -1.4176, -1.2207, -1.3887, -1.4276, -1.0926, -1.2531, -1.3538, -1.2761, -1.1706, -1.3323, -1.3072, -1.3627, -1.0738, -1.2761, -1.2971, -1.1243, -1.2662, -1.0263, -1.1336, -0.6919, -1.0611, -1.0452, -0.9676, -0.8375, -1.1052, -0.8617, -0.6976, -0.9645, -0.8217, -0.7416, -0.7680, -0.7487, -0.6404, -0.6585, -0.5160, -0.5734, -0.4544, -0.3548, -0.5550, -0.2872, -0.5568, -0.5904, -0.3976, -0.4900, -0.1741, -0.2722, -0.2806, -0.1337, -0.0358, -0.0823, -0.0529, -0.0329, -0.0633, -0.0557, -0.1239, -0.1643, 0.0480, 0.1552, -0.0435, 0.2747, 0.0024, 0.2003, 0.0367, 0.3316, 0.1624, 0.2403, 0.3370, 0.1532, 0.3668, 0.2304, 0.4972, 0.5929, 0.3871, 0.6119, 0.5865, 0.6421, 0.5096, 0.6150, 0.5723, 0.6167, 0.7797, 0.6231, 0.6377, 0.9428, 0.5556, 0.7645, 0.6638, 0.4265, 0.4847, 0.8218, 0.3660, 0.4582, 0.5175, -0.6659, -0.4388, -0.8430, -0.7341, -0.4414, -0.5791, -0.6557, -0.7276, -0.4589, -0.6036, -0.4564, -0.6255, -0.6331, -0.6575, -0.4523, -0.5745, -0.4909, -0.4819, -0.4261, -0.4568, -0.3690, -0.4626, -0.5536, -0.2953, -0.5536, -0.5103, -0.2853, -0.2833, -0.1254, -0.3125, -0.0550, -0.3311, -0.1740, -0.1884, -0.1742, -0.0724, 0.0436, 0.0017, -0.0240, 0.2086, 0.1771, -0.1345, 0.0630, 0.1153, 0.2607, 0.3638, 0.2768, 0.2395, 0.4097, 0.3963, 0.4683, 0.5509, 0.4934, 0.4535, 0.5158, 0.3567, 0.5306, 0.6044, 0.6605, 0.6521, 0.6950, 0.7203, 0.6228, 0.6046, 0.8399, 0.8635, 0.9783, 0.8162, 0.8786, 0.8792, 1.2167, 1.0865, 1.1895, 1.0693, 1.1609, 0.9985, 1.0642, 0.9084, 1.2205, 1.1302, 1.2652, 1.1474, 0.9671, 1.2496, 1.1915, 1.2055, 1.2922, 1.5113, 1.3225, 1.4352, 1.5476, 1.1648, 1.3783, 1.4360, 1.4545, 1.2762, 1.3074, 1.3631, 1.3704, 1.5025, /*Second row*/ -0.2140, -0.2108, -0.1714, -0.0356, 0.2183, 0.2104, -0.0088, -0.0112, 0.2611, 0.3050, 0.3084, 0.2285, 0.3889, 0.3172, 0.2294, 0.3511, 0.2285, 0.1800, 0.3830, 0.4024, 0.4343, 0.4546, 0.5531, 0.4366, 0.5352, 0.6903, 0.5706, 0.6211, 0.3925, 0.6801, 0.5915, 0.6738, 0.7937, 0.4094, 0.9375, 0.8461, 0.7222, 0.6609, 0.6901, 0.9520, 0.8394, 0.7916, 0.7570, 0.5260, 1.0404, 0.8415, 0.9159, 0.9733, 0.9597, 0.7515, 0.7913, 0.7874, 0.9387, 0.8282, 0.6716, 0.8793, 0.7465, 0.7691, 0.6191, 0.9164, 0.8220, 0.7640, 0.6195, 0.6749, 0.8239, 0.5973, 0.7304, 0.6996, 0.5767, 0.6471, 0.4028, 0.6390, 0.9071, 0.5980, 0.5627, 0.5919, 0.6838, 0.4179, 0.5818, 0.4410, 0.3200, 0.4471, 0.4527, 0.3645, 0.1701, 0.3644, 0.1532, 0.4038, 0.2053, 0.0350, -0.0498, 0.1661, -0.0729, 0.0675, -0.3096, -0.1284, 0.1117, -0.1848, -0.0372, 0.0718, -0.0879, -0.1076, 0.0384, 0.0066, -0.0214, -0.1831, -0.1936, -0.0980, 0.0592, -0.0522, -0.0738, 0.0059, -0.0071, -0.2125, -0.4056, -0.3016, -0.3823, -0.2313, -0.4488, -0.7583, -0.5339, -0.5700, -0.4820, -0.4901, -0.5232, -0.5142, -0.4817, -0.8765, -0.5628, -0.7358, -0.6668, -0.7392, -0.7508, -0.7142, -0.9030, -0.8873, -0.7846, -0.8245, -0.6875, -0.7454, -0.7149, -0.7254, -0.7262, -0.6192, -0.7144, -0.9667, -0.6916, -0.9017, -0.9633, -0.8407, -0.8531, -0.8954, -0.9079, -0.8683, -0.9853, -0.8038, -0.9842, -0.9084, -0.6880, -0.7589, -0.6349, -0.8472, -0.9383, -0.9943, -0.7895, -0.8038, -0.5556, -0.6824, -0.7625, -0.7662, -0.5972, -0.5120, -0.7345, -0.5620, -0.3713, -0.6260, -0.6160, -0.5891, -0.4883, -0.3904, -0.3449, -0.4375, -0.2620, -0.0903, -0.4001, -0.0510, -0.2635, -0.2463, -0.2388, -0.0223, -0.0940, -0.1269, -0.0965, -0.0159, 0.0792, 0.0848, -0.1011, 0.0876, 0.3170, 0.1342 };
//static double v2[400] = { /*First row*/ -1.3583, -1.5049, -1.5156, -1.1801, -1.3588, -1.2245, -1.0528, -1.3870, -1.3284, -1.0855, -1.5378, -1.2338, -1.3296, -1.1645, -1.3227, -1.2013, -1.3009, -1.1555, -1.2047, -1.1679, -0.8691, -0.9678, -1.1191, -1.0332, -1.1343, -1.0792, -1.1052, -0.9497, -0.9198, -0.9520, -0.9452, -0.9207, -0.8876, -0.6510, -0.5037, -0.8894, -0.6996, -0.4753, -0.8198, -0.5498, -0.5836, -0.5970, -0.4536, -0.7281, -0.6171, -0.4407, -0.4368, -0.3296, -0.6192, -0.5108, -0.5266, -0.5112, -0.2432, -0.0858, -0.4430, -0.3339, -0.1368, -0.1030, 0.0664, 0.1782, -0.0518, -0.0634, 0.0181, 0.1302, 0.0648, 0.1191, 0.0597, 0.0611, 0.1422, 0.2678, 0.2814, 0.1818, 0.4043, 0.1779, 0.3511, 0.1773, 0.3872, 0.6803, 0.3249, 0.4425, 0.5761, 0.4672, 0.4377, 0.6880, 0.5569, 0.8432, 0.7974, 0.4596, 0.6612, 0.7845, 0.6533, 0.6062, 0.5570, 0.4803, 0.5068, 0.6994, 0.5847, 0.5731, 0.6385, 0.5697, -0.7624, -0.5350, -0.7049, -0.7145, -0.7448, -0.7579, -0.6022, -0.6539, -0.6093, -0.5181, -0.3276, -0.7968, -0.7073, -0.2268, -0.2542, -0.6220, -0.4609, -0.6222, -0.3290, -0.4211, -0.3797, -0.4087, -0.3476, -0.3610, -0.3479, -0.1508, -0.5372, -0.4187, -0.0789, -0.2672, -0.2082, -0.3269, 0.1315, 0.0042, -0.1201, -0.0762, -0.3049, 0.0733, 0.0423, -0.0422, 0.1095, 0.1851, -0.1527, 0.1558, 0.2106, 0.3293, 0.2611, 0.1503, 0.2482, 0.3656, 0.3096, 0.4578, 0.5571, 0.6960, 0.6284, 0.5579, 0.3849, 0.8484, 0.6694, 0.9665, 1.0687, 0.5584, 0.9050, 0.9331, 0.6966, 0.8271, 0.8938, 0.7596, 0.9966, 0.7985, 0.8517, 0.9981, 1.1322, 1.0717, 0.8644, 0.9947, 1.1702, 1.1638, 0.9658, 1.3382, 1.0576, 1.1382, 1.2054, 1.6196, 1.2550, 1.2911, 1.2076, 1.3668, 1.2338, 1.4051, 1.2932, 1.2042, 1.2087, 1.4852, 1.3221, 1.5053, 1.5650, 1.3629, 1.3534, 1.4212, /*Second row*/ -0.1474, 0.0836, -0.1002, 0.0641, -0.2777, -0.1711, 0.0988, 0.1492, -0.0113, 0.2532, 0.0717, 0.1682, 0.2915, 0.2440, 0.1394, 0.5083, 0.6296, 0.4168, 0.1416, 0.5305, 0.5327, 0.4387, 0.3135, 0.4191, 0.4624, 0.5064, 0.4872, 0.5650, 0.5185, 0.3563, 0.4836, 0.5023, 0.8515, 0.7416, 0.7863, 0.9049, 0.8588, 0.7928, 0.7809, 0.7086, 0.8581, 0.7989, 0.7171, 0.6368, 0.9799, 0.9114, 0.9009, 0.8352, 0.8213, 0.9554, 0.9092, 0.7130, 0.9075, 1.0810, 0.8878, 0.7926, 0.7463, 0.9377, 1.0203, 0.7238, 0.8931, 0.7257, 0.6909, 0.7976, 0.6252, 0.7312, 0.7420, 0.3709, 0.6202, 0.7639, 0.6382, 0.6423, 0.8075, 0.6428, 0.7160, 0.6323, 0.5238, 0.3594, 0.3179, 0.6203, 0.3667, 0.3347, 0.4224, 0.3487, 0.2565, 0.4796, 0.0242, 0.1543, 0.0785, 0.2180, 0.3741, 0.3275, 0.0384, -0.0178, 0.1574, -0.2140, -0.0821, -0.1030, -0.1444, -0.3283, 0.1250, 0.0550, 0.1439, -0.2224, -0.1390, 0.0773, 0.0067, 0.1490, -0.0313, -0.1131, 0.0179, -0.2429, -0.1800, -0.1726, -0.4180, -0.1482, -0.3423, -0.4155, -0.3951, -0.2809, -0.3341, -0.4015, -0.4883, -0.6606, -0.4843, -0.6760, -0.4048, -0.5894, -0.5518, -0.5558, -0.7294, -0.5967, -0.6921, -0.8267, -0.8548, -0.7411, -0.9591, -0.7843, -0.7887, -0.7038, -0.8230, -0.8315, -0.8418, -0.7263, -0.6371, -0.7524, -0.7220, -0.7463, -0.9144, -0.7840, -0.7365, -0.7480, -0.6960, -0.7344, -0.7179, -0.8860, -0.9248, -0.7322, -0.8129, -0.7213, -0.6632, -0.7540, -0.7249, -0.8975, -0.6129, -0.8301, -0.7588, -0.5104, -0.6345, -0.8534, -0.4594, -0.6043, -0.6231, -0.4942, -0.5582, -0.3807, -0.4713, -0.4109, -0.6387, -0.4284, -0.5362, -0.4388, -0.3789, -0.4907, -0.3726, -0.0927, -0.1597, -0.0939, -0.1217, -0.2198, -0.0717, 0.1372, -0.3007, -0.0667, -0.1304, 0.0538, 0.0199, 0.2726, 0.0620, 0.2353 };

static double v1[400] = { /*First row*/ 1.5025, -1.2887, -1.3371, -1.3821, -1.2948, -1.5536, -1.5409, -1.5559, -1.4176, -1.2207, -1.3887, -1.4276, -1.0926, -1.2531, -1.3538, -1.2761, -1.1706, -1.3323, -1.3072, -1.3627, -1.0738, -1.2761, -1.2971, -1.1243, -1.2662, -1.0263, -1.1336, -0.6919, -1.0611, -1.0452, -0.9676, -0.8375, -1.1052, -0.8617, -0.6976, -0.9645, -0.8217, -0.7416, -0.7680, -0.7487, -0.6404, -0.6585, -0.5160, -0.5734, -0.4544, -0.3548, -0.5550, -0.2872, -0.5568, -0.5904, -0.3976, -0.4900, -0.1741, -0.2722, -0.2806, -0.1337, -0.0358, -0.0823, -0.0529, -0.0329, -0.0633, -0.0557, -0.1239, -0.1643, 0.0480, 0.1552, -0.0435, 0.2747, 0.0024, 0.2003, 0.0367, 0.3316, 0.1624, 0.2403, 0.3370, 0.1532, 0.3668, 0.2304, 0.4972, 0.5929, 0.3871, 0.6119, 0.5865, 0.6421, 0.5096, 0.6150, 0.5723, 0.6167, 0.7797, 0.6231, 0.6377, 0.9428, 0.5556, 0.7645, 0.6638, 0.4265, 0.4847, 0.8218, 0.3660, 0.4582, 0.5175, -0.6659, -0.4388, -0.8430, -0.7341, -0.4414, -0.5791, -0.6557, -0.7276, -0.4589, -0.6036, -0.4564, -0.6255, -0.6331, -0.6575, -0.4523, -0.5745, -0.4909, -0.4819, -0.4261, -0.4568, -0.3690, -0.4626, -0.5536, -0.2953, -0.5536, -0.5103, -0.2853, -0.2833, -0.1254, -0.3125, -0.0550, -0.3311, -0.1740, -0.1884, -0.1742, -0.0724, 0.0436, 0.0017, -0.0240, 0.2086, 0.1771, -0.1345, 0.0630, 0.1153, 0.2607, 0.3638, 0.2768, 0.2395, 0.4097, 0.3963, 0.4683, 0.5509, 0.4934, 0.4535, 0.5158, 0.3567, 0.5306, 0.6044, 0.6605, 0.6521, 0.6950, 0.7203, 0.6228, 0.6046, 0.8399, 0.8635, 0.9783, 0.8162, 0.8786, 0.8792, 1.2167, 1.0865, 1.1895, 1.0693, 1.1609, 0.9985, 1.0642, 0.9084, 1.2205, 1.1302, 1.2652, 1.1474, 0.9671, 1.2496, 1.1915, 1.2055, 1.2922, 1.5113, 1.3225, 1.4352, 1.5476, 1.1648, 1.3783, 1.4360, 1.4545, 1.2762, 1.3074, 1.3631, 1.3704, /*Second row*/ 0.1342, -0.2140, -0.2108, -0.1714, -0.0356, 0.2183, 0.2104, -0.0088, -0.0112, 0.2611, 0.3050, 0.3084, 0.2285, 0.3889, 0.3172, 0.2294, 0.3511, 0.2285, 0.1800, 0.3830, 0.4024, 0.4343, 0.4546, 0.5531, 0.4366, 0.5352, 0.6903, 0.5706, 0.6211, 0.3925, 0.6801, 0.5915, 0.6738, 0.7937, 0.4094, 0.9375, 0.8461, 0.7222, 0.6609, 0.6901, 0.9520, 0.8394, 0.7916, 0.7570, 0.5260, 1.0404, 0.8415, 0.9159, 0.9733, 0.9597, 0.7515, 0.7913, 0.7874, 0.9387, 0.8282, 0.6716, 0.8793, 0.7465, 0.7691, 0.6191, 0.9164, 0.8220, 0.7640, 0.6195, 0.6749, 0.8239, 0.5973, 0.7304, 0.6996, 0.5767, 0.6471, 0.4028, 0.6390, 0.9071, 0.5980, 0.5627, 0.5919, 0.6838, 0.4179, 0.5818, 0.4410, 0.3200, 0.4471, 0.4527, 0.3645, 0.1701, 0.3644, 0.1532, 0.4038, 0.2053, 0.0350, -0.0498, 0.1661, -0.0729, 0.0675, -0.3096, -0.1284, 0.1117, -0.1848, -0.0372, 0.0718, -0.0879, -0.1076, 0.0384, 0.0066, -0.0214, -0.1831, -0.1936, -0.0980, 0.0592, -0.0522, -0.0738, 0.0059, -0.0071, -0.2125, -0.4056, -0.3016, -0.3823, -0.2313, -0.4488, -0.7583, -0.5339, -0.5700, -0.4820, -0.4901, -0.5232, -0.5142, -0.4817, -0.8765, -0.5628, -0.7358, -0.6668, -0.7392, -0.7508, -0.7142, -0.9030, -0.8873, -0.7846, -0.8245, -0.6875, -0.7454, -0.7149, -0.7254, -0.7262, -0.6192, -0.7144, -0.9667, -0.6916, -0.9017, -0.9633, -0.8407, -0.8531, -0.8954, -0.9079, -0.8683, -0.9853, -0.8038, -0.9842, -0.9084, -0.6880, -0.7589, -0.6349, -0.8472, -0.9383, -0.9943, -0.7895, -0.8038, -0.5556, -0.6824, -0.7625, -0.7662, -0.5972, -0.5120, -0.7345, -0.5620, -0.3713, -0.6260, -0.6160, -0.5891, -0.4883, -0.3904, -0.3449, -0.4375, -0.2620, -0.0903, -0.4001, -0.0510, -0.2635, -0.2463, -0.2388, -0.0223, -0.0940, -0.1269, -0.0965, -0.0159, 0.0792, 0.0848, -0.1011, 0.0876, 0.3170 };
static double v2[400] = { /*First row*/ 1.4212, -1.3583, -1.5049, -1.5156, -1.1801, -1.3588, -1.2245, -1.0528, -1.3870, -1.3284, -1.0855, -1.5378, -1.2338, -1.3296, -1.1645, -1.3227, -1.2013, -1.3009, -1.1555, -1.2047, -1.1679, -0.8691, -0.9678, -1.1191, -1.0332, -1.1343, -1.0792, -1.1052, -0.9497, -0.9198, -0.9520, -0.9452, -0.9207, -0.8876, -0.6510, -0.5037, -0.8894, -0.6996, -0.4753, -0.8198, -0.5498, -0.5836, -0.5970, -0.4536, -0.7281, -0.6171, -0.4407, -0.4368, -0.3296, -0.6192, -0.5108, -0.5266, -0.5112, -0.2432, -0.0858, -0.4430, -0.3339, -0.1368, -0.1030, 0.0664, 0.1782, -0.0518, -0.0634, 0.0181, 0.1302, 0.0648, 0.1191, 0.0597, 0.0611, 0.1422, 0.2678, 0.2814, 0.1818, 0.4043, 0.1779, 0.3511, 0.1773, 0.3872, 0.6803, 0.3249, 0.4425, 0.5761, 0.4672, 0.4377, 0.6880, 0.5569, 0.8432, 0.7974, 0.4596, 0.6612, 0.7845, 0.6533, 0.6062, 0.5570, 0.4803, 0.5068, 0.6994, 0.5847, 0.5731, 0.6385, 0.5697, -0.7624, -0.5350, -0.7049, -0.7145, -0.7448, -0.7579, -0.6022, -0.6539, -0.6093, -0.5181, -0.3276, -0.7968, -0.7073, -0.2268, -0.2542, -0.6220, -0.4609, -0.6222, -0.3290, -0.4211, -0.3797, -0.4087, -0.3476, -0.3610, -0.3479, -0.1508, -0.5372, -0.4187, -0.0789, -0.2672, -0.2082, -0.3269, 0.1315, 0.0042, -0.1201, -0.0762, -0.3049, 0.0733, 0.0423, -0.0422, 0.1095, 0.1851, -0.1527, 0.1558, 0.2106, 0.3293, 0.2611, 0.1503, 0.2482, 0.3656, 0.3096, 0.4578, 0.5571, 0.6960, 0.6284, 0.5579, 0.3849, 0.8484, 0.6694, 0.9665, 1.0687, 0.5584, 0.9050, 0.9331, 0.6966, 0.8271, 0.8938, 0.7596, 0.9966, 0.7985, 0.8517, 0.9981, 1.1322, 1.0717, 0.8644, 0.9947, 1.1702, 1.1638, 0.9658, 1.3382, 1.0576, 1.1382, 1.2054, 1.6196, 1.2550, 1.2911, 1.2076, 1.3668, 1.2338, 1.4051, 1.2932, 1.2042, 1.2087, 1.4852, 1.3221, 1.5053, 1.5650, 1.3629, 1.3534, /*Second row*/ 0.2353, -0.1474, 0.0836, -0.1002, 0.0641, -0.2777, -0.1711, 0.0988, 0.1492, -0.0113, 0.2532, 0.0717, 0.1682, 0.2915, 0.2440, 0.1394, 0.5083, 0.6296, 0.4168, 0.1416, 0.5305, 0.5327, 0.4387, 0.3135, 0.4191, 0.4624, 0.5064, 0.4872, 0.5650, 0.5185, 0.3563, 0.4836, 0.5023, 0.8515, 0.7416, 0.7863, 0.9049, 0.8588, 0.7928, 0.7809, 0.7086, 0.8581, 0.7989, 0.7171, 0.6368, 0.9799, 0.9114, 0.9009, 0.8352, 0.8213, 0.9554, 0.9092, 0.7130, 0.9075, 1.0810, 0.8878, 0.7926, 0.7463, 0.9377, 1.0203, 0.7238, 0.8931, 0.7257, 0.6909, 0.7976, 0.6252, 0.7312, 0.7420, 0.3709, 0.6202, 0.7639, 0.6382, 0.6423, 0.8075, 0.6428, 0.7160, 0.6323, 0.5238, 0.3594, 0.3179, 0.6203, 0.3667, 0.3347, 0.4224, 0.3487, 0.2565, 0.4796, 0.0242, 0.1543, 0.0785, 0.2180, 0.3741, 0.3275, 0.0384, -0.0178, 0.1574, -0.2140, -0.0821, -0.1030, -0.1444, -0.3283, 0.1250, 0.0550, 0.1439, -0.2224, -0.1390, 0.0773, 0.0067, 0.1490, -0.0313, -0.1131, 0.0179, -0.2429, -0.1800, -0.1726, -0.4180, -0.1482, -0.3423, -0.4155, -0.3951, -0.2809, -0.3341, -0.4015, -0.4883, -0.6606, -0.4843, -0.6760, -0.4048, -0.5894, -0.5518, -0.5558, -0.7294, -0.5967, -0.6921, -0.8267, -0.8548, -0.7411, -0.9591, -0.7843, -0.7887, -0.7038, -0.8230, -0.8315, -0.8418, -0.7263, -0.6371, -0.7524, -0.7220, -0.7463, -0.9144, -0.7840, -0.7365, -0.7480, -0.6960, -0.7344, -0.7179, -0.8860, -0.9248, -0.7322, -0.8129, -0.7213, -0.6632, -0.7540, -0.7249, -0.8975, -0.6129, -0.8301, -0.7588, -0.5104, -0.6345, -0.8534, -0.4594, -0.6043, -0.6231, -0.4942, -0.5582, -0.3807, -0.4713, -0.4109, -0.6387, -0.4284, -0.5362, -0.4388, -0.3789, -0.4907, -0.3726, -0.0927, -0.1597, -0.0939, -0.1217, -0.2198, -0.0717, 0.1372, -0.3007, -0.0667, -0.1304, 0.0538, 0.0199, 0.2726, 0.0620 };

static matrix data[2] = { {.data = v1, .w = 200, .h = 2}, {.data = v2, .w = 200, .h = 2} };

int main(int argc, char *argv[]) {
    getchar();
    
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    gmc_result r = gmc(data, 2, 2, 1.0d, 0);

    gettimeofday(&end, 0);
    if(r.cluster_num != 2) printf("Couldn't find requested cluster number (%d). Got %d clusters\n", 2, r.cluster_num);
    printf("Iteration %d: λ=%lf\n", r.iterations, r.lambda);
    printf("Time measured: %.3f seconds.\n", end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6);
    
    //for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
    //printf("\n");
    //print(r.U);

    return 0;
}
