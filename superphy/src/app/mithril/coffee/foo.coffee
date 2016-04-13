# coffeelint: disable=max_line_length

class FooFoo
    @controller: (args) ->
        args.id = args.id || 0
        @not_saved = m.prop("")
        @saved = m.prop(localStorage.getItem("saved#{args.id}"))
        @onunload = () ->
            localStorage.setItem("saved#{args.id}", @saved())
        return @

    @view: (ctrl) ->
        m('.'
            m('hr')
            
            "Saved"
            m('input[type=text]', {
                oninput: m.withAttr('value', ctrl.saved)
                value: ctrl.saved()
            })
            "Not Saved"
            m('input[type=text]', {
                oninput: m.withAttr('value', ctrl.not_saved)
                value: ctrl.not_saved()
            })
        )