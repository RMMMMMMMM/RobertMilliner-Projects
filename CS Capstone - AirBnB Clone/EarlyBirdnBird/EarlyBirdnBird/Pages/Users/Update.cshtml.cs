using Infrastructure.Models;
using Infrastructure.Interfaces;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;

namespace EarlyBirdnBird.Pages.Users
{
    public class UpdateModel : PageModel
    {
        private readonly IUnitOfWork _unitOfWork;
        private readonly UserManager<AppUser> _userManager;
        private readonly RoleManager<IdentityRole> _roleManager;

        public UpdateModel(IUnitOfWork unitOfWork, UserManager<AppUser> userManager, RoleManager<IdentityRole> roleManager)
        {
            _unitOfWork = unitOfWork;
            _userManager = userManager;
            _roleManager = roleManager;

        }

        [BindProperty]
        public AppUser AppUser { get; set; }
        public List<string> UserRoles { get; set; }
        public List<string> AllRoles { get; set; }
        public List<string> OldRoles { get; set; }

        public async Task OnGetAsync(string id)
        {
            AppUser = _unitOfWork.AppUser.Get(u => u.Id == id);
            var roles = await _userManager.GetRolesAsync(AppUser);
            UserRoles = roles.ToList();
            OldRoles = roles.ToList();
            AllRoles = _roleManager.Roles.Select(r => r.Name).ToList();

        }

        public async Task<IActionResult> OnPostAsync()
        {
            var newRoles = Request.Form["roles"];
            UserRoles = newRoles.ToList();
            var OldRoles = await _userManager.GetRolesAsync(AppUser); //ones in DB
            var rolesToAdd = new List<string>();
            var user = _unitOfWork.AppUser.Get(u => u.Id == AppUser.Id);

            user.FName = AppUser.FName;
            user.LName = AppUser.LName;
            user.Email = AppUser.Email;
            user.PhoneNumber = AppUser.PhoneNumber;
            _unitOfWork.AppUser.Update(user);
            _unitOfWork.Commit();

            //update their roles
            foreach (var r in UserRoles)
            {
                if (!OldRoles.Contains(r)) // new Role
                {
                    rolesToAdd.Add(r);
                }
            }

            foreach(var r in OldRoles)
            {
                if(!UserRoles.Contains(r)) // remove
                {
                    var result = await _userManager.RemoveFromRoleAsync(user, r);
                }
            }

            var result1 = await _userManager.AddToRolesAsync(user, rolesToAdd.AsEnumerable());
            return RedirectToPage("./Index", new { success = true, message = "Update Successful" });
        }
    }
}
