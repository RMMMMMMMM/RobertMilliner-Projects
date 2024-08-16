using Infrastructure.Interfaces;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore.Migrations.Internal;
using Infrastructure.Models;
using DataAccess;
using Microsoft.EntityFrameworkCore;

namespace EarlyBirdnBird.Pages.Users
{
    public class IndexModel : PageModel
    {
        private readonly IUnitOfWork _unitOfWork;
        private readonly UserManager<AppUser> _userManager;
        private readonly ApplicationDbContext _dbContext;

        public IndexModel(IUnitOfWork unitOfWork, UserManager<AppUser> userManager, ApplicationDbContext db)
        {
            _userManager = userManager;
            _unitOfWork = unitOfWork;
            _dbContext = db;
        }

        public IEnumerable<AppUser> ApplicationUsers { get; set; }
        public Dictionary<string, List<string>> UserRoles { get; set; }
        public bool Success {  get; set; }
        public string Message { get; set; }

        //for searching with email
        [BindProperty(SupportsGet = true)]
        public string EmailSearch { get; set; }

        //for filtering by role
        [BindProperty(SupportsGet = true)]
        public string RoleFilter { get; set; }

        //pagination properties
        [BindProperty(SupportsGet = true)]
        public int CurrentPage { get; set; } = 1;
        private int PageSize { get; } = 20;
        public int TotalPages { get; set; }

        public async Task OnGetAsync(bool success = false, string message = null)
        {
            Success = success;
            Message = message;
            UserRoles = new Dictionary<string, List<string>>();

            var usersQuery = _userManager.Users.AsQueryable();

            if (!string.IsNullOrEmpty(EmailSearch))
            {
                usersQuery = usersQuery.Where(x => x.Email.Contains(EmailSearch));
            }

            var usersList = await usersQuery.ToListAsync();
            foreach (var user in usersList)
            {
                var userRoles = await _userManager.GetRolesAsync(user);
                UserRoles[user.Id] = userRoles.ToList();
            }

            if (!string.IsNullOrEmpty(RoleFilter))
            {
                usersList = usersList.Where(user => UserRoles[user.Id].Contains(RoleFilter)).ToList();
            }

            var userCount = usersList.Count;
            TotalPages = (int)Math.Ceiling(userCount / (double)PageSize);

            CurrentPage = Math.Max(CurrentPage, 1);
            CurrentPage = Math.Min(CurrentPage, TotalPages);

            ApplicationUsers = usersList.Skip((CurrentPage - 1) * PageSize).Take(PageSize).ToList();
        }



        public async Task<IActionResult> OnPostLockUnlock(string id)
        {
            var user = _unitOfWork.AppUser.Get(u => u.Id == id);
            if(user.LockoutEnd == null) //unlocked
            {
                user.LockoutEnd = DateTime.Now.AddYears(100);
                user.LockoutEnabled = true;
            }
            else if(user.LockoutEnd > DateTime.Now)//unlock
            {
                user.LockoutEnd = DateTime.Now;
                user.LockoutEnabled = false;
            }
            else
            {
                user.LockoutEnd = DateTime.Now.AddYears(100);
                user.LockoutEnabled = true;
            }
            _unitOfWork.AppUser.Update(user);
            await _unitOfWork.CommitAsync();
            return RedirectToPage();
        }

        public async Task<IActionResult> OnPostDeleteAsync(string id)
        {
            var user = await _userManager.FindByIdAsync(id);
            if (user != null)
            {
                var result = await _userManager.DeleteAsync(user);
                if (result.Succeeded)
                {
                    Success = true;
                    Message = "User has been deleted successfully.";
                }
                else
                {
                    Success = false;
                    Message = "Error deleting user.";
                }
            }
            else
            {
                Success = false;
                Message = "User not found.";
            }
            return RedirectToPage(new { success = Success, message = Message });
        }

    }
}
