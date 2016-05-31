# coffeelint: disable=max_line_length

dragdrop = (element, options) ->
    activate = (e) ->
        e.preventDefault()
        return
    deactivate = ->
    update = (e) ->
        e.preventDefault()
        if typeof options.onchange == 'function'
            options.onchange (e.dataTransfer or e.target).files
        return
    options = options or {}
    element.addEventListener 'dragover', activate
    element.addEventListener 'dragleave', deactivate
    element.addEventListener 'dragend', deactivate
    element.addEventListener 'drop', deactivate
    element.addEventListener 'drop', update
    window.addEventListener 'blur', deactivate
    return

submodule = (module, args) ->
    module.view.bind this, new (module.controller)(args)

class Foo
    @upload: (files) ->
        formData = new FormData
        i = 0
        while i < files.length
            formData.append('file' + i, files[i])
            i++
        response = m.request(
            method: 'POST'
            url: "http://#{location.hostname}:5000/data/"
            data: formData
            datatype:"multipart/form-data"
            serialize: (value) ->
                value
            type: (response) ->
            console.log("RESPONSE: ", JSON.stringify(response))
        )
        console.log("Formdata: " + (JSON.stringify(formData.entries())))

    @controller: (options) ->
        options = options or {}
        @onchange = options.onchange
        return @
    @view = (ctrl) ->
        m('.uploader'
            {

            }
            "FOO"
        )

class Uploader extends Page
    Routes.add('/uploader', @)
    Routes.add('/upload', @)

    @controller: ->
        @title= -> 'Upload something'
        @foo = m.prop()
        return @
    @view: (ctrl) ->
        m('.'
            m('h1', ctrl.title())
            m('input[type=file]', {
                oninput:m.withAttr('value', ctrl.foo)
                value:ctrl.foo()
                placeholder:"Username"
            })
            m('input[type=button]'
                value:"PUSH ME"
                onchange: (files) ->
                    ctrl.files = files
                onclick: (files) ->
                    alert(JSON.stringify(ctrl.foo))
            )
        )



###
        m('.'
            m('h1', ctrl.title())
            m.component(Uploader, onchange: (files) ->
                Uploader.upload(files)
                #.then(->
                #    alert 'uploaded!'
                #    return
                #)
                return
            )
            m('p', 'drop a file in the box above')
        )
###