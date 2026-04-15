/* Hide edit this page in apidocs */
window.onload = function(){ removeLinks(); };

function removeLinks() {
    var githubPath = document.querySelector("div.tocsection.editthispage a").href;
    var dirNames = githubPath.split('/');
    if (dirNames.includes('_autosummary') == true ) {
        var tocSection = document.querySelector("div.tocsection.editthispage");
        tocSection.style.display = "none";
        document.querySelector("li.nav-item > a[data-bs-original-title='Launch on Binder']").closest("li").style.display = "none"
    }
}