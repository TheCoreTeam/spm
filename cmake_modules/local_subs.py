"""
 @file local_subs.py

 Python SPM specific substitution rules for the Precision Generator script.

 @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.2.3
 @author Mathieu Faverge
 @date 2023-12-11

"""
subs = {
    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'normal': [
        # pattern                single                  double                  single-complex          double-complex
        #'12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890')
        ('int',                  'float',                'double',               'spm_complex32_t',     r'\bspm_complex64_t'   ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmComplex32',        r'\bSpmComplex64'      ),
        ('SpmPattern',           'SpmFloat',             'SpmDouble',            'SpmFloat',            r'\bSpmDouble'         ),
        ('MPI_INT',              'MPI_FLOAT',            'MPI_DOUBLE',           'MPI_COMPLEX32',        'MPI_COMPLEX64'       ),

        # ----- Variables
        (r'\b',                 r'szero\b',             r'dzero\b',             r'czero\b',             r'zzero\b'             ),
        (r'\b',                 r'sone\b',              r'done\b',              r'cone\b',              r'zone\b'              ),

        # ----- SPM Prefixes
        ('spm_p',                'spm_s',                'spm_d',                'spm_c',                'spm_z'               ),
    ]
}

exceptfrom = []
