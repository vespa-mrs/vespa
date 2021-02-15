
"""\
FloatSpinMultiplier widget configuration

@copyright: 2016-2019 Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

# In the 'events' section below, we have deliberately added a preceding
# space. This is a hack to avoid the inherent wcodegen.PythonMixin.cn()
# method that checks the first 4 chars for 'EVT_' and then prefixes
# a 'wx.' in front. For this widget, we only want the actual event 
# alias from the floatspin_multiplier_base module, EVT_FLOATSPIN. The prefix
# space keeps wxglade from decorating it.
#
#    'events': {
#        ' EVT_FLOATSPIN': {
#            'type': 'FloatSpin',
#        },
#    },
#

# keep synchronous between wxSpinCtrl and wxTextCtrl
config = {
    'wxklass': 'FloatSpinMultiplier',
    'style_defs': {
        'wxSP_ARROW_KEYS': { 'desc': _('The user can use arrow keys to change the value.') },
        'wxSP_WRAP':       { 'desc': _('The value wraps at the minimum and maximum.') },
        'wxALIGN_LEFT': {
            'desc': _('Same as wxTE_LEFT for wxTextCtrl: the text is left aligned.'),
            'exclude': 'wxALIGN_CENTRE_HORIZONTAL|wxALIGN_RIGHT',
            'supported_by': ('wx3',) },
        'wxALIGN_CENTRE_HORIZONTAL': {
            'desc': _('Same as wxTE_CENTRE for wxTextCtrl: the text is centered.'),
            'exclude': 'wxALIGN_LEFT|wxALIGN_RIGHT',
            'supported_by': ('wx3',) },
        'wxALIGN_RIGHT': {
            'desc': _('Same as wxTE_RIGHT for wxTextCtrl: the text is right aligned (this is the default).'),
            'exclude': 'wxALIGN_LEFT|wxALIGN_CENTRE_HORIZONTAL',
            'supported_by': ('wx3',) },
        'wxTE_PROCESS_ENTER': {
            'desc': _('The control will generate the event wxEVT_TEXT_ENTER (otherwise pressing Enter key is '
                      'either processed internally by the control or used for navigation between dialog controls).'), },
        'wxTE_PROCESS_TAB':  {'obsolete':True},
        'wxTE_MULTILINE':    {'obsolete':True},
        'wxTE_PASSWORD':     {'obsolete':True},
        'wxTE_READONLY':     {'obsolete':True},
        'wxTE_RICH':         {'obsolete':True},
        'wxTE_RICH2':        {'obsolete':True},
        'wxTE_AUTO_URL':     {'obsolete':True},
        'wxTE_NOHIDESEL':    {'obsolete':True},
        'wxHSCROLL':         {'obsolete':True},
        'wxTE_NO_VSCROLL':   {'obsolete':True},
        'wxTE_LEFT':         {'rename_to':'wxALIGN_LEFT'},
        'wxTE_CENTRE':       {'rename_to':'wxALIGN_CENTRE_HORIZONTAL'},
        'wxTE_RIGHT':        {'rename_to':'wxALIGN_RIGHT'},
        'wxTE_DONTWRAP':     {'obsolete':True},
        'wxTE_LINEWRAP':     {'obsolete':True},
        'wxTE_CHARWRAP':     {'obsolete':True},
        'wxTE_WORDWRAP':     {'obsolete':True},
        'wxTE_BESTWRAP':     {'obsolete':True},
        'wxTE_CAPITALIZE':   {'obsolete':True}
    },
    'default_style': 'wxSP_ARROW_KEYS',
    'style_list': ['wxSP_ARROW_KEYS', 'wxSP_WRAP',
                   'wxTE_PROCESS_ENTER',
                   'wxALIGN_LEFT', 'wxALIGN_CENTRE_HORIZONTAL', 'wxALIGN_RIGHT'],
    'events': {
        ' EVT_FLOATSPIN': {'type': 'wxFloatSpin',},
    },
}
