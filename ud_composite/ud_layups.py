# -*- coding: utf-8 -*-
# =============================================================
#  ud_layups.py — Layup catalog (74 candidates)
#
#  Source: ud_ply_angle_candidates.md
#
#  Symmetric notation [a/b]s expanded to explicit [a, b, b, a].
#
#  Each entry:
#    {
#       'id'    : int,         # layup ID (1..74)
#       'group' : '2-ply' | '3-ply' | '4-ply' | '6-ply' | '8-ply',
#       'label' : str,         # original notation, e.g. '[45/-45]s'
#       'angles': list of int, # explicit ply angles, outermost-first
#       'notes' : str,
#    }
#
#  Ply thickness = 0.23 mm (CU250).
# =============================================================

PLY_THICKNESS = 0.23  # mm

# ----- 2-ply (t = 0.46 mm) -----
_2PLY = [
    ( 1, '[0/90]'   , [0, 90]   , 'cross-ply'),
    ( 2, '[0/0]'    , [0, 0]    , 'unidirectional axial'),
    ( 3, '[90/90]'  , [90, 90]  , 'unidirectional transverse'),
    ( 4, '[10/-10]' , [10, -10] , 'angle-ply very narrow'),
    ( 5, '[15/-15]' , [15, -15] , 'angle-ply narrow'),
    ( 6, '[20/-20]' , [20, -20] , 'angle-ply'),
    ( 7, '[25/-25]' , [25, -25] , 'angle-ply'),
    ( 8, '[30/-30]' , [30, -30] , 'angle-ply mid'),
    ( 9, '[35/-35]' , [35, -35] , 'angle-ply'),
    (10, '[40/-40]' , [40, -40] , 'angle-ply'),
    (11, '[45/-45]' , [45, -45] , 'woven 45 equivalent'),
    (12, '[50/-50]' , [50, -50] , 'angle-ply'),
    (13, '[55/-55]' , [55, -55] , 'angle-ply'),
    (14, '[60/-60]' , [60, -60] , 'angle-ply wide'),
    (15, '[65/-65]' , [65, -65] , 'angle-ply'),
    (16, '[70/-70]' , [70, -70] , 'angle-ply'),
    (17, '[75/-75]' , [75, -75] , 'angle-ply near transverse'),
    (18, '[0/45]'   , [0, 45]   , 'asymmetric'),
    (19, '[0/60]'   , [0, 60]   , 'asymmetric'),
    (20, '[0/30]'   , [0, 30]   , 'asymmetric'),
]

# ----- 3-ply (t = 0.69 mm) -----
_3PLY = [
    (21, '[0/90/0]'    , [0, 90, 0]    , 'cross-ply axial outer'),
    (22, '[90/0/90]'   , [90, 0, 90]   , 'cross-ply transverse outer'),
    (23, '[45/0/-45]'  , [45, 0, -45]  , 'quasi-iso variant'),
    (24, '[45/90/-45]' , [45, 90, -45] , '45 series + transverse'),
    (25, '[60/0/-60]'  , [60, 0, -60]  , 'wide + axial'),
    (26, '[30/0/-30]'  , [30, 0, -30]  , 'narrow + axial'),
    (27, '[0/45/0]'    , [0, 45, 0]    , 'axial dominant + shear'),
    (28, '[0/60/0]'    , [0, 60, 0]    , 'axial dominant + wide'),
    (29, '[45/0/45]'   , [45, 0, 45]   , 'asymmetric shear'),
    (30, '[30/60/-30]' , [30, 60, -30] , 'mixed angle'),
    (31, '[0/45/90]'   , [0, 45, 90]   , 'spread asymmetric'),
    (32, '[20/0/-20]'  , [20, 0, -20]  , 'narrow + axial'),
    (33, '[40/0/-40]'  , [40, 0, -40]  , 'mid + axial'),
    (34, '[50/0/-50]'  , [50, 0, -50]  , 'wide + axial'),
    (35, '[70/0/-70]'  , [70, 0, -70]  , 'very wide + axial'),
]

# ----- 4-ply (t = 0.92 mm) -----
_4PLY = [
    (36, '[45/-45]s'        , [45, -45, -45, 45]  , 'woven 2-ply equivalent, symmetric'),
    (37, '[45/-45/45/-45]'  , [45, -45, 45, -45]  , 'repeated'),
    (38, '[0/90]s'          , [0, 90, 90, 0]      , 'cross-ply symmetric'),
    (39, '[30/-30]s'        , [30, -30, -30, 30]  , 'mid symmetric'),
    (40, '[60/-60]s'        , [60, -60, -60, 60]  , 'wide symmetric'),
    (41, '[20/-20]s'        , [20, -20, -20, 20]  , 'narrow symmetric'),
    (42, '[40/-40]s'        , [40, -40, -40, 40]  , 'mid symmetric'),
    (43, '[50/-50]s'        , [50, -50, -50, 50]  , 'mid symmetric'),
    (44, '[70/-70]s'        , [70, -70, -70, 70]  , 'wide symmetric'),
    (45, '[75/-75]s'        , [75, -75, -75, 75]  , 'near transverse symmetric'),
    (46, '[0/45/-45/0]'     , [0, 45, -45, 0]     , 'axial outer + shear inner'),
    (47, '[45/0/0/-45]'     , [45, 0, 0, -45]     , 'shear outer + axial inner'),
    (48, '[0/90/90/0]'      , [0, 90, 90, 0]      , 'cross-ply symmetric duplicate'),
    (49, '[90/0/0/90]'      , [90, 0, 0, 90]      , 'transverse outer'),
    (50, '[0/60/-60/0]'     , [0, 60, -60, 0]     , 'axial outer + wide'),
    (51, '[0/30/-30/0]'     , [0, 30, -30, 0]     , 'axial outer + narrow'),
    (52, '[30/60/-60/-30]'  , [30, 60, -60, -30]  , 'mixed angle'),
    (53, '[45/90/-45/0]'    , [45, 90, -45, 0]    , 'quasi-iso asymmetric'),
    (54, '[0/45/90/-45]'    , [0, 45, 90, -45]    , 'quasi-iso asymmetric 2'),
]

# ----- 6-ply (t = 1.38 mm) -----
_6PLY = [
    (55, '[45/0/-45]s'   , [45, 0, -45, -45, 0, 45]  , 'quasi-iso symmetric'),
    (56, '[60/0/-60]s'   , [60, 0, -60, -60, 0, 60]  , 'wide symmetric'),
    (57, '[30/0/-30]s'   , [30, 0, -30, -30, 0, 30]  , 'narrow symmetric'),
    (58, '[0/45/90]s'    , [0, 45, 90, 90, 45, 0]    , 'full spread symmetric'),
    (59, '[0/45/-45]s'   , [0, 45, -45, -45, 45, 0]  , 'axial + shear symmetric'),
    (60, '[90/45/-45]s'  , [90, 45, -45, -45, 45, 90], 'transverse + shear symmetric'),
    (61, '[0/60/-60]s'   , [0, 60, -60, -60, 60, 0]  , 'axial + wide symmetric'),
    (62, '[0/30/-30]s'   , [0, 30, -30, -30, 30, 0]  , 'axial + narrow symmetric'),
    (63, '[45/-45/0]s'   , [45, -45, 0, 0, -45, 45]  , 'shear outer symmetric'),
    (64, '[45/-45/90]s'  , [45, -45, 90, 90, -45, 45], 'shear + transverse symmetric'),
    (65, '[20/0/-20]s'   , [20, 0, -20, -20, 0, 20]  , 'very narrow symmetric'),
    (66, '[40/0/-40]s'   , [40, 0, -40, -40, 0, 40]  , 'mid symmetric'),
    (67, '[50/0/-50]s'   , [50, 0, -50, -50, 0, 50]  , 'mid-wide symmetric'),
    (68, '[70/0/-70]s'   , [70, 0, -70, -70, 0, 70]  , 'very wide symmetric'),
]

# ----- 8-ply (t = 1.84 mm) -----
_8PLY = [
    (69, '[45/0/-45/90]s'   , [45, 0, -45, 90, 90, -45, 0, 45]   , 'quasi-iso standard'),
    (70, '[45/-45/0/90]s'   , [45, -45, 0, 90, 90, 0, -45, 45]   , 'quasi-iso variant'),
    (71, '[0/45/90/-45]s'   , [0, 45, 90, -45, -45, 90, 45, 0]   , 'quasi-iso variant 2'),
    (72, '[30/60/-60/-30]s' , [30, 60, -60, -30, -30, -60, 60, 30], 'mixed angle symmetric'),
    (73, '[0/45/-45/0]s'    , [0, 45, -45, 0, 0, -45, 45, 0]     , 'axial outer symmetric'),
    (74, '[45/0/0/-45]s'    , [45, 0, 0, -45, -45, 0, 0, 45]     , 'shear outer symmetric'),
]


def _pack(group_name, raw):
    out = []
    for lid, label, angles, notes in raw:
        out.append({
            'id'    : lid,
            'group' : group_name,
            'label' : label,
            'angles': list(angles),
            'notes' : notes,
        })
    return out


LAYUPS = (
    _pack('2-ply', _2PLY) +
    _pack('3-ply', _3PLY) +
    _pack('4-ply', _4PLY) +
    _pack('6-ply', _6PLY) +
    _pack('8-ply', _8PLY)
)


def get_layup(layup_id):
    """Look up a layup record by ID (1..74)."""
    for lay in LAYUPS:
        if lay['id'] == layup_id:
            return lay
    raise KeyError('Unknown layup_id: %s' % layup_id)


def layup_thickness(layup_id):
    return len(get_layup(layup_id)['angles']) * PLY_THICKNESS


if __name__ == '__main__':
    print('Total layups: %d' % len(LAYUPS))
    by_group = {}
    for lay in LAYUPS:
        by_group.setdefault(lay['group'], []).append(lay['id'])
    for g in ('2-ply', '3-ply', '4-ply', '6-ply', '8-ply'):
        ids = by_group.get(g, [])
        print('  %s : %d entries (IDs %s..%s, t=%.2f mm)' % (
            g, len(ids), ids[0] if ids else '-',
            ids[-1] if ids else '-',
            (len(get_layup(ids[0])['angles']) * PLY_THICKNESS) if ids else 0))
