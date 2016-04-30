"""
class C1 extends ComponentTemplate
    @controller: (args) ->
        value = m.prop("")
        view: () ->
            return m("div"
                m("button", {onclick: args.arg}, "Custom Message")
                m("input[type=text]", {
                    oninput: m.withAttr('value', value), value: value()
                })
                m("button", {onclick: -> alert(value())}, "Submit")
            )

class C2 extends C1
    @view: (ctrl) ->
        return m('b', ctrl.view())

class BoldC1 extends ComponentTemplate
    @controller: (args) ->
        view: () ->
            m('b', m.component(C1, {arg: args.arg}))

class Input extends ComponentTemplate
    @controller: (args) ->
        view: () ->
            return m("input[type=text]", {
                oninput: m.withAttr('value', args.value), value: args.value()
            })

class GroupBrowse2 extends PageTemplate
    @controller: (args) ->
        value = m.prop("")
        view: () ->
            return m(".page",
                super()
                m.component(Input, {value: value})
                value()
            )
"""

class Input extends ComponentTemplate
    @controller: (args) ->
        view: () ->
            return m("input[type=text]", {
                oninput: m.withAttr('value', args.value), value: args.value()
            })
