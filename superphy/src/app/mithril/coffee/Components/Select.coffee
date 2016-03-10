class Select2
    view: (ctrl, attrs) ->
        selectedId = attrs.value().id
        m('select', {config: Select2.config(attrs)}, [ attrs.data.map((item) ->
            args = value: item.id
            #    Set selected option
            if item.id == selectedId
                args.selected = 'selected'
            return m('option', args, item.name)
        )])

    config: (ctrl) ->
        return (element, isInitialized) -> 
            if typeof jQuery is not 'undefined' and typeof jQuery.fn.select2 is not 'undefined'
                el = $(element)
                if not isInitialized
                    el.select2().on 'change', (e) ->
                        id = el.select2('val')
                        m.startComputation()

                        # Set value to selected option
                        ctrl.data.map (d) ->
                            if d.id = id then ctrl.value(d)
                            return
                        if typeof ctrl.onchange == 'function'
                            ctrl.onchance(el.select2('val'))

                        m.endComputation()
                        return
                el.val(ctrl.value().id).trigger('change')
            else
                console.warn('ERROR: You need jquery and Select22 in the page')





